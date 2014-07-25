import pickle
import spectrum_utilities

result =  spectrum_utilities.compute_sla_power_spectrum(
    spectrum_utilities.regions, dset='unfiltered')

pickle.dump(result, open('global_AVISO_data.pkl','w'))