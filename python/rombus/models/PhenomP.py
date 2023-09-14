from typing import NamedTuple

import lal  # type: ignore
import lalsimulation  # type: ignore
import numpy as np

from rombus.model import RombusModel


class Model(RombusModel):

    
    # Set some constants
    l_1 = 0
    l_2 = 0
    f_min = 20.0
    f_max = 1024.0
    delta_F = 1.0 / 4.0
    n_f = int((f_max - f_min) / delta_F) + 1
    f_min_index = int(f_min / delta_F)
    WFdict = lal.CreateDict()

    # N.B.: mypy struggles with NamedTuples, so typing is turned off for the following
    params.add("m1", 30, 35)  # type: ignore # noqa F821
    params.add("m2", 30, 35)  # type: ignore # noqa F821
    params.add("chi1L", 0, 0.1)  # type: ignore # noqa F821
    params.add("chi2L", 0, 0.1)  # type: ignore # noqa F821
    params.add("chip", 0, 0.1)  # type: ignore # noqa F821
    params.add("thetaJ", -np.pi / 2, np.pi / 2)  # type: ignore # noqa F821
    params.add("alpha", 0, np.pi / 2)  # type: ignore # noqa F821
    params.set_validation(lambda p: p.m1 >= p.m2)  # type: ignore # noqa F821

    def m1m2_to_Mc(m1,m2):
    """Chirp mass from m1, m2"""
    return (m1*m2)**(3./5.)/(m1+m2)**(1./5.)

    def earth_motion_time_delay(geocent_time,chirp_mass,frequency):
    '''
    Eq. just below (9) in https://arxiv.org/pdf/1710.05325.pdf,
    note that the equation in the paper above assumes G=c=1
    '''
        t_f = geocent_time - (5/256)*(chirp_mass*solar_mass)**(-5/3) * \
            (np.pi*frequency)**(-8/3)*speed_of_light**5 * gravitational_constant**(-5/3)
        return t_f

    
    tc = None
    Mc_min = m1m2_to_Mc(self.model_params['m1'].min, self.model_params['m2'].min)

    times = [earth_motion_time_delay(tc, Mc, 2**(i+1))- earth_motion_time_delay(tc, Mc, fmin)]

    fmin_ind = int(f_min%2)
    fmax_ind = int(f_max%2)
    
    for i in range(fmin_ind, fmax_ind+1):
        times.append(
            earth_motion_time_delay(tc, Mc, 2**(i+1)) - earth_motion_time_delay(tc, Mc, 2**(i)
                                                                               )
    delta_fs = 1/times

    freq_range = np.linspace(fmin, 2**(fmin_ind+1), int((2**(fmin_ind+1)-fmin)/delta_f[0])+1)

    for i in range(fmin_ind+1, fmax_ind)+1):
        freq_range = np.append(freq_range, np.linspace(2**i, 2**(i+1), int((2**i)/delta_f[i-fmin_ind+1)])))
    
    # Set the domain over-and-on which the ROM will be defined
    coordinate.set("f", f_min, f_max, n_f, label="$f$", dtype=np.dtype("float64"), values=freq_range)  # type: ignore # noqa F821

    # Set the ordinate the model will map the domain to
    ordinate.set("h", label="$h$", dtype=complex)  # type: ignore # noqa F821
    

    def compute(self, params: NamedTuple, domain: np.ndarray) -> np.ndarray:

        lalsimulation.SimInspiralWaveformParamsInsertTidalLambda1(self.WFdict, self.l_1)
        lalsimulation.SimInspiralWaveformParamsInsertTidalLambda2(self.WFdict, self.l_2)
        if not np.array_equiv(domain, self.domain):
            h = lalsimulation.SimIMRPhenomPFrequencySequence(
                domain,
                params.chi1L,  # type: ignore
                params.chi2L,  # type: ignore
                params.chip,  # type: ignore
                params.thetaJ,  # type: ignore
                params.m1 * lal.lal.MSUN_SI,  # type: ignore
                params.m2 * lal.lal.MSUN_SI,  # type: ignore
                1e6 * lal.lal.PC_SI * 100,
                params.alpha,  # type: ignore
                0,
                40,
                lalsimulation.IMRPhenomPv2NRTidal_V,
                lalsimulation.NRTidalv2_V,
                self.WFdict,
            )
            h = h[0].data.data
        else:
            h = lalsimulation.SimIMRPhenomP(
                params.chi1L,  # type: ignore
                params.chi2L,  # type: ignore
                params.chip,  # type: ignore
                params.thetaJ,  # type: ignore
                params.m1 * lal.lal.MSUN_SI,  # type: ignore
                params.m2 * lal.lal.MSUN_SI,  # type: ignore
                1e6 * lal.lal.PC_SI * 100,
                params.alpha,  # type: ignore
                0,
                self.delta_F,
                self.f_min,
                self.f_max,
                40,
                lalsimulation.IMRPhenomPv2NRTidal_V,
                lalsimulation.NRTidalv2_V,
                self.WFdict,
            )
            h = h[0].data.data[self.f_min_index : len(h[0].data.data)]
            if len(h) < self.n_domain:
                h = np.append(h, np.zeros(self.n_domain - len(h), dtype=complex))

        return h
