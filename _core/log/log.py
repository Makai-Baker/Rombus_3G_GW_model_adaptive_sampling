"""This module provides a `LogStream` class for generating logging
information.  It is intended for the generation of course-grained reporting of
program execution for the user and should not be used in performance-critical
situations, in tight loops, etc.

Formatting is organized by indenting levels which can be
increased/decreased by calling the open/close methods of the stream
respectively.
"""

# For legacy-Python compatibility
from __future__ import print_function
from functools import update_wrapper
import traceback

from typing import List

import os
import sys
import time
import datetime

import inspect
from functools import wraps

intervals = (
    ("weeks", 604800),  # 60 * 60 * 24 * 7
    ("days", 86400),  # 60 * 60 * 24
    ("hours", 3600),  # 60 * 60
    ("minutes", 60),
    ("seconds", 1),
)

string_types = (str, bytes)

TIMING_AUTO_SECONDS_DEFAULT = 0.5


def is_nonstring_iterable(object_in):
    """Determine if an object is a non-string iterable.

    :param object: An object of any type.
    :return: Boolean.  True if object is a non-string iterable.
    """
    return hasattr(object_in, "__iter__") and not isinstance(object_in, string_types)


def format_time(seconds, granularity=None):
    """Create a nice ASCII representation of a time interval, given in seconds.

    :param seconds: Time in seconds
    :param granularity: The maximum number of interval levels to report
    :return: string
    """
    result = []

    if seconds < 1:
        if seconds < 1e-3:
            result.append(f"{int(1e6*seconds)} microseconds")
        else:
            result.append(f"{int(1e3*seconds)} milliseconds")
    else:
        for i_interval, [name, count] in enumerate(intervals):
            value = seconds // count
            if value:
                seconds -= value * count
                if name == intervals[-1][0] or i_interval == granularity:
                    result.append("%.1f %s" % (value, name))
                else:
                    if value == 1:
                        name = name.rstrip("s")
                    result.append("%d %s" % (value, name))

    if granularity:
        result = ", ".join(result[:granularity])
    else:
        result = ", ".join(result)

    # Replace the last ',' with ' and'
    result_split = result.rsplit(",", 1)
    if len(result_split) > 1:
        result = " and".join(result_split)

    return result


class LogStream(object):
    """This class provides a file pointer for logging user feedback and methods
    for writing to it."""

    def __init__(
        self,
        fp_out=None,
        verbosity=True,
        n_indent_max=10,
        exception_handler=None,
        time_elapsed_auto_seconds=TIMING_AUTO_SECONDS_DEFAULT,
    ):
        """
        :param fp_out: An optional file pointer to use for the log.
        :param verbosity: An optional parameter that sets the default verbosity of the stream.
        :param n_indent_max: maximum number of logging levels to keep track of.  Anything exceeding this is not printed.
        """
        # File pointer where the stream will write to
        self.set_fp(fp_out)

        # Number of spaces to indent for each indent-level
        self.indent_size = 3

        # Set the maximum number of indent levels to render
        self.n_indent_max = n_indent_max

        if exception_handler is None:
            self._exception_handler = self.handle_exception
        else:
            self._exception_handler = exception_handler

        self.time_elapsed_auto_seconds = time_elapsed_auto_seconds

        # These lists will have one entry per indent-level
        self.t_last = [time.time()]
        self.n_lines = [0]
        self.splice = [None]
        self.time_elapsed = ["auto"]

        # This list will be a stack with one entry per verbosity state.  Initialize with the given default.
        self.verbosity: List = []
        self.verbosity_default = verbosity
        self.set_verbosity(self.verbosity_default)

        # Indicates whether the last-written line
        # ended with a new line
        self.hanging = False

        self._halt = False

    def open(self, msg, splice=None, time_elapsed="auto"):
        """Open a new indent bracket for the log.

        :param msg: An object with a __str__ method, or a list thereof
        :return: None
        """
        self._print(msg + "...", unhang=True, indent=True)
        self.t_last.append(time.time())
        self.n_lines.append(0)
        self.splice.append(splice)
        self.time_elapsed.append(time_elapsed)
        if splice:
            self._splice_line(splice, True)

    def close(self, msg=None):
        """Close a new indent bracket for the log.

        :param msg: An object with a __str__ method, or a list thereof
        :return: None
        """

        # Sanity checks
        if self._n_indent() < 1:
            self.error(Exception("Invalid log closure."))

        # Decrement the indent level and fetch the info about the level we are closing
        t_last = self.t_last.pop()
        n_lines = self.n_lines.pop()
        splice = self.splice.pop()
        time_elapsed = self.time_elapsed.pop()

        # This must be called every time because we need the
        # pop on t_last to keep track of the indenting level
        dt = time.time() - t_last

        if time_elapsed == "auto":
            if dt >= self.time_elapsed_auto_seconds:
                time_elapsed = True
            else:
                time_elapsed = False

        if splice:
            self._splice_line(splice, False)

        # Generate message
        if msg is not None:
            msg_time = ""
            if time_elapsed:
                dt_txt = format_time(dt)
                if len(dt_txt) > 0:
                    msg_time = " (%s)" % (dt_txt)
            self._print(msg + msg_time, unhang=(n_lines > 1))
        self._unhang()

    def handle_exception(self, exception: Exception) -> None:

        # Make sure to unhang the logger
        self._unhang()
        self.blankline()

        # Print exception and stack trace then exit
        traceback.print_exception(exception)

        # Halt the logger, in case we're inside nested log.contexts
        self.halt()

        raise

    def callable(
        self,
        msg=None,
        dump_args=False,
        dump_returns=False,
        time_elapsed="auto",
        default_verbosity="unset",
    ):
        """Decorator to add in-bound and out-bound logging to a callable.

        :param msg: string describing in-bound logging message
        :param dump_args: log calling arguments if True
        :param dump_returns: log returned values if True
        :param time_elapsed: log ellapsed time if True
        :param default_verbosity: default verbosity
        :return: decorated callable
        """

        def decorated_callable(func):
            if default_verbosity != "unset":
                func = self.add_verbosity(default=default_verbosity)(func)

            @wraps(func)
            def wrapper(*args, **kwargs):

                # Print function call
                if msg:
                    self.open(msg, time_elapsed=time_elapsed)
                else:
                    self.open(
                        "Calling %s.%s()..." % (func.__module__, func.__qualname__),
                        time_elapsed=time_elapsed,
                    )

                # Report arguments
                if dump_args:
                    self.open("Inputs:")
                    func_args = inspect.signature(func).bind(*args, **kwargs).arguments
                    if len(func_args) > 1:
                        for i, item in enumerate(func_args.items()):
                            self.comment("{} = {!r}".format(*item))
                    else:
                        self.comment("None")
                    self.close()

                # Run method
                if dump_args or dump_returns:
                    self.open("Running...")
                r = func(*args, **kwargs)
                if dump_args or dump_returns:
                    self.close("Done.")

                # Report returns
                if dump_returns and r is not None:
                    if len(r) > 1:
                        for i, r_i in enumerate(r):
                            self.comment("Return[%d]: %s" % (i, r_i))
                    else:
                        self.comment("Return: " + r)

                # We don't need to display time elapsed twice ...
                self.close("Done.")
                return r

            return update_wrapper(wrapper, func)

        return decorated_callable

    class _Context:
        def __init__(self, stream, *args, **kwargs):
            self.stream = stream
            self.args = args
            self.kwargs = kwargs
            self._exception_handler = stream._exception_handler

        def __enter__(self):
            self.stream.open(*self.args, **self.kwargs)

        def __exit__(self, exc_type, exc_val, exc_tb):
            if exc_type is not None:
                self._exception_handler(exc_val)
            else:
                self.stream.close("Done.")
            return True

    def context(self, *args, **kwargs):
        return LogStream._Context(self, *args, **kwargs)

    class _Progress:
        def __init__(self, stream, msg, target, *args, reverse=False, **kwargs):
            self.stream = stream
            self.msg = msg
            self.target = target
            self.args = args
            self.kwargs = kwargs
            self._exception_handler = stream._exception_handler
            self.reverse = reverse

            self.i_update = 1
            self.n_update = 10
            if self.reverse:
                self.report_next = self.target * (self.n_update - 1) / self.n_update
            else:
                self.report_next = self.target / self.n_update

        def reset_next(self, progress):
            if self.reverse:
                self.report_next = max(
                    [
                        self.target,
                        progress
                        - (progress - self.target) / (self.n_update - self.i_update),
                    ]
                )
            else:
                self.report_next = min(
                    [
                        self.target,
                        progress
                        + (self.target - progress) / (self.n_update - self.i_update),
                    ]
                )

        def update(self, progress):
            if self.i_update < self.n_update:
                if self.reverse:
                    if progress <= self.report_next:
                        progress_pc = int(min([100, 100 * progress / self.target]))
                        self.stream.comment(f"{progress_pc:2}% complete.")
                        self.report_next = max(
                            [
                                self.target,
                                progress
                                - (progress - self.target)
                                / (self.n_update - self.i_update),
                            ]
                        )
                        self.i_update = self.i_update + 1
                else:
                    if progress >= self.report_next:
                        progress_pc = int(min([100, 100 * progress / self.target]))
                        self.stream.comment(f"{progress_pc:2}% complete.")
                        self.report_next = min(
                            [
                                self.target,
                                progress
                                + (self.target - progress)
                                / (self.n_update - self.i_update),
                            ]
                        )
                        self.i_update = self.i_update + 1

        def __enter__(self):
            self.stream.open(self.msg, *self.args, **self.kwargs)
            return self

        def __exit__(self, exc_type, exc_val, exc_tb):
            if exc_type is not None:
                self._exception_handler(exc_val)
            else:
                self.stream.close("Done.")
            return True

    def progress(self, *args, **kwargs):
        return LogStream._Progress(self, *args, **kwargs)

    def test(
        self,
        msg=None,
        dump_args=False,
        dump_returns=False,
        time_elapsed=False,
        default_verbosity="unset",
    ):
        """Decorator to add in-bound and out-bound logging to a test.

        :param msg: string describing in-bound logging message
        :param dump_args: log calling arguments if True
        :param dump_returns: log returned values if True
        :param time_elapsed: log ellapsed time if True
        :param default_verbosity: default verbosity
        :return: decorated test
        """
        self.__init__(fp_out=sys.stdout)
        self.blankline()
        return self.callable(
            msg=msg,
            dump_args=dump_args,
            dump_returns=dump_returns,
            time_elapsed=time_elapsed,
            default_verbosity=default_verbosity,
        )

    def methods(
        self,
        dump_args=True,
        dump_returns=True,
        time_elapsed=True,
        default_verbosity="unset",
    ):
        """Decorator for automating the logging of all method calls of a class.

        :param dump_args: log calling arguments if True
        :param dump_returns: log returned values if True
        :param time_elapsed: log ellapsed time if True
        :param default_verbosity: default verbosity
        :return: decorated method
        """

        def decorated_class_declaration(Cls):
            class NewCls(object):
                def __init__(new_self, *args, **kwargs):
                    new_self.oInstance = Cls(*args, **kwargs)

                def __getattribute__(new_self, s):
                    """this is called whenever any attribute of a NewCls object
                    is accessed.

                    This function first tries to get the attribute off
                    NewCls. If it fails then it tries to fetch the
                    attribute from self.oInstance (an instance of the
                    decorated class). If it manages to fetch the
                    attribute from self.oInstance, and the attribute is
                    an instance method then the method decorator is
                    applied.
                    """
                    try:
                        x = super(NewCls, new_self).__getattribute__(s)
                    except AttributeError:
                        pass
                    else:
                        return x
                    x = new_self.oInstance.__getattribute__(s)
                    if isinstance(x, type(new_self.__init__)):
                        return self.callable(
                            dump_args=dump_args,
                            dump_returns=dump_returns,
                            time_elapsed=time_elapsed,
                            default_verbosity=default_verbosity,
                        )(x)
                    else:
                        return x

            return NewCls

        return decorated_class_declaration

    def set_fp(self, fp_out=None):
        """Set the file pointer to be used for logging.  Default is `sys.stdout`.

        :param fp_out: File pointer
        :return: None
        """

        if "PYTEST_CURRENT_TEST" in os.environ:
            self.fp = open(os.devnull, "w")
        else:
            # self.fp = sys.stderr
            self.fp = sys.stdout

    def set_verbosity(self, verbosity=True):
        """Add a new (and make it current) verbosity state to the stream's
        stack of verbosity states.

        This method takes either a boolean flag indicating whether logging is active, or an integer indicating
        the maximum indenting level that will be rendered.  It can be removed using the
        :py:meth:`~.log.LogStream.unset_verbosity` method.  See the
        :py:meth:`~.log.LogStream.check_verbosity` method for an account of how the verbosity passed to this
        method is interpreted.

        :param verbosity: A boolean flag indicating if logging is active, or an integer indicating the verbosity level
        :return: None
        """

        # Check validity of the given verbosity
        if type(verbosity) not in (bool, int):
            self.error(
                TypeError(
                    "Invalid datatype {%s} being added to log stream's verbosity state."
                    % (type(verbosity))
                )
            )

        if not isinstance(verbosity, bool):
            verbosity += self._n_indent() - 1

        # Add a state to the stack
        self.verbosity.append(verbosity)

    def unset_verbosity(self):
        """Revert stream to a previous verbosity state if one exists; the
        default state otherwise.

        :return: None
        """
        if len(self.verbosity) > 0:
            self.verbosity.pop()

    def verbosity_level(self, verbosity):
        """Convert a verbosity state value to a corresponding verbosity level.

        :param verbosity: Verbosity state value
        :return: Integer indent level
        """

        # Default result
        result = self.n_indent_max

        # If state is a bool and evaluates to false, return -1 (i.e. a value always > self._n_indent()
        if isinstance(verbosity, bool):
            if not verbosity:
                result = -1

        # ... else, if it's an integer, return it or n_indent_max
        elif isinstance(verbosity, int):
            result = min([verbosity, self.n_indent_max])

        # ... else, unsupported data type ... throw an error
        else:
            self.error(
                Exception(
                    "Can not interpret verbosity level of a verbosity state with unsupported type {%s}."
                    % (type(verbosity))
                )
            )

        return result

    def check_verbosity(self):
        """Check if the stream is active.

        :return: A boolean indicating if rendering is active on the stream
        """
        # If the verbosity stack is empty, use the default
        if len(self.verbosity) < 1:
            max_active_level = self.verbosity_level(self.verbosity_default)
        else:
            max_active_level = self.n_indent_max
            for state in self.verbosity:
                max_active_level = min([max_active_level, self.verbosity_level(state)])

        return max_active_level >= self._n_indent()

    def add_verbosity(self, default=True):
        """Decorator to add a 'verbosity' parameter - and the functions to implement it - to a callable.

        :param default: default verbosity
        :return: decorated callable
        """

        def decorated_callable(func):
            if default == "unset":
                return func
            else:

                @wraps(func)
                def wrapper(*args, **kwargs):
                    # Fetch callable's verbosity; emit a TypeError if not present
                    if "verbosity" not in kwargs:
                        verbosity = default
                    else:
                        # Check the signature to make sure that we aren't passing
                        # verbosity to a function that does not declare it and
                        # does not support kwargs
                        (
                            func_args,
                            func_varargs,
                            func_keywords,
                            func_defaults,
                        ) = inspect.getargspec(func)
                        if "verbosity" not in func_args and not func_keywords:
                            verbosity = kwargs.pop("verbosity")
                        else:
                            verbosity = kwargs.get("verbosity")

                    # Set callable's verbosity
                    self.set_verbosity(verbosity)

                    # Run callable
                    r = func(*args, **kwargs)

                    # Unset callable's verbosity
                    self.unset_verbosity()

                    return r

                return update_wrapper(wrapper, func)

        return decorated_callable

    def comment(
        self,
        msg,
        unhang=True,
        overwrite=False,
        blankline_before=False,
        blankline_after=False,
    ):
        """Add a one-line comment to the log.

        :param msg: An object with a __str__ method, or a list thereof
        :param unhang:
        :param overwrite:
        :return: None
        """
        if blankline_before:
            self.blankline()
        self._print(msg, unhang=unhang, indent=True, overwrite=overwrite)
        if blankline_after:
            self.blankline()

    def append(self, msg):
        """Add to the end of the current line in the log.

        :param msg: An object with a __str__ method, or a list thereof
        :return: None
        """
        self._print(msg + "...", unhang=False, indent=False)

    def progress_bar(self, gen, count, *args, **kwargs):
        """Display a progress bar for a generator.

        :param gen: Generator
        :param count: Number of generator iterations
        :param args: Positional arguments to pass to the generator
        :param kwargs: Keyword arguments to pass to the generator
        :return: None
        """

        # Initialize counter
        width = 30
        msg_len_last = 0
        start_time = time.time()
        self.comment("[%s] Remaining:" % (" " * width), unhang=True)

        # Iterate
        for iteration, result in enumerate(gen(*args, **kwargs)):
            fraction_complete = float(iteration + 1) / float(count)
            ticks = int(fraction_complete * float(width + 1))
            secs_elapsed = time.time() - start_time
            secs_estimate = int(secs_elapsed / fraction_complete)
            secs_remaining = secs_estimate - secs_elapsed
            if secs_remaining > 0:
                msg = "[%s%s] Remaining: %s" % (
                    "#" * ticks,
                    " " * (width - ticks),
                    str(datetime.timedelta(seconds=secs_remaining)).split(".")[0],
                )
                msg_len = len(msg)

                # Make sure to blank-out any old underlying text
                if msg_len < msg_len_last:
                    msg += " " * (msg_len_last - msg_len)
                msg_len_last = msg_len
                self.comment(msg, unhang=False, overwrite=True)

        # Finalize counter
        msg = "[%s%s] Time elapsed: %s" % (
            "#" * ticks,
            " " * (width - ticks),
            str(datetime.timedelta(seconds=secs_elapsed)).split(".")[0],
        )
        msg_len = len(msg)
        if msg_len < msg_len_last:
            msg += " " * (msg_len_last - msg_len)
        self.comment(msg, unhang=False, overwrite=True)

    def blankline(self):
        """Print a blank line to the stream.

        :return: None
        """
        self.comment("\n", unhang=True)

    def raw(self, msg):
        """Print raw, unformatted text to the log.

        :param msg: An object with a __str__ method, or a list thereof
        :return: None
        """
        self._print(msg, unhang=True, indent=False)

    def halt(self):
        self._unhang()
        self._halt = True

    def _splice_line(self, splice_msg, flag_start):
        """Create splice lines in the log for isolating sections of the stream.

        This method is intended to be used when uncontrolled output from other sources are polluting the stream.  Open an indentation
        block around cases like this using the splice keyword argument, and a clearly identifiable line will be
        rendered at the start and end of the section.

        :param splice_msg:
        :param flag_start:
        :return:
        """
        n_splice = 40
        n_lead_min = 10
        lead_char = "="
        msg = " " + splice_msg + " - "
        if flag_start:
            msg += "start "
        else:
            msg += "end "
        n_msg = len(msg)
        n_lead = int((n_splice - len(msg)) / 2)
        if n_lead <= 0:
            n_splice = n_msg + 2 * n_lead_min
            n_lead = n_lead_min
            n_tail = n_lead_min
        else:
            n_tail = n_splice - n_msg - n_lead
        self._print(
            n_lead * lead_char + msg + n_tail * lead_char + "\n",
            unhang=True,
            indent=False,
        )

    def _print_to_fp(self, msg, **kwargs):

        if not self._halt:
            print(msg, **kwargs)
            self.fp.flush()

    def _print(
        self,
        msg,
        unhang=True,
        indent=True,
        overwrite=False,
        iterables_allowed=True,
        **kwargs,
    ):
        """This method is the main driver of output to the stream, but should
        be accessed through other methods.

        :param msg: An object with a __str__ method, or a list thereof
        :param unhang: Boolean flag indicating whether to start with a carriage return
        :param indent: Boolean flag indicating whether to start the line with an indent
        :param iterables_allowed: Boolean flag indicating whether to accept an iterable msg
        :param kwargs: keyword arguments to be passed to the print function
        :return: None
        """
        # Check if rendering is active on the stream
        if self.check_verbosity():

            # Optionally unhang the stream
            if unhang:
                self._unhang()

            # This will fail for strings but pass for lists, etc.
            if is_nonstring_iterable(msg):
                if overwrite:
                    self.error(
                        Exception("Log stream overwriting not permitted for iterables.")
                    )
                if not iterables_allowed:
                    self.error(
                        Exception(
                            "An iterable was passed to a log stream method which does not accept them."
                        )
                    )
                for line in msg:
                    self._print(line, indent=indent, overwrite=overwrite, **kwargs)
            # ... render a non-iterable object ...
            else:
                # If msg is a string (or converts to one) with newline characters, break-it-up
                # and recall this method with the result to treat it as an iterable
                str_msg = str(msg)
                msg_split = str_msg.splitlines(True)
                if len(msg_split) > 1:
                    self._print(
                        msg_split,
                        indent=indent,
                        overwrite=overwrite,
                        iterables_allowed=iterables_allowed,
                        **kwargs,
                    )
                # ... else, render a single line
                else:
                    if not self.hanging and len(str_msg) > 0:
                        self.n_lines[-1] += 1
                    if overwrite or (not self.hanging and indent):
                        self._indent(overwrite=overwrite)
                    self._print_to_fp(str_msg, end="", file=self.fp, **kwargs)
                    if str_msg.endswith("\n"):
                        self.hanging = False
                    else:
                        self.hanging = True

    def _unhang(self):
        """If the log did not previously end with a newline, add one.

        :return: None
        """
        if self.hanging:
            self._print_to_fp("")
            self.n_lines[-1] += 1
            self.hanging = False

    def _indent(self, overwrite=False):
        """Write the appropriate indent for this line (with an option to
        overwrite)

        :param overwrite: Boolean flag indicating whether to overwrite the current line
        :return: None
        """
        if overwrite:
            self._print_to_fp("\r", end="", file=self.fp)
        self._print_to_fp(
            self.indent_size * self._n_indent() * " ", end="", file=self.fp
        )

    def _n_indent(self):
        """Return the current indent level of the stream.

        :return: Integer
        """
        return len(self.t_last) - 1
