import pytest


def pytest_collection_modifyitems(config, items):
    keywordexpr = config.option.keyword
    markexpr = config.option.markexpr
    if keywordexpr or markexpr:
        return  # let pytest handle this

    skip_lalsuite = pytest.mark.skip(reason="lalsuite not selected")
    for item in items:
        if "lalsuite" in item.keywords:
            item.add_marker(skip_lalsuite)
