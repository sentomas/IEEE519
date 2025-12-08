module IEEE519API

using ..IEEE519
using CSV, DataFrames

"""
Analyze a CSV file for IEEE519 compliance.
Arguments:
- csv_filename: path to CSV file
- fundamental_freq: 50 or 60 Hz
- pcc_kv: PCC voltage in kV
- plot_results: whether to plot results (default true)
Returns: NamedTuple with voltage, current, sampling_freq, nyquist_ok, compliance
"""
function analyze_csv(csv_filename::String; fundamental_freq::Float64=50.0, pcc_kv::Float64=12.47, plot_results::Bool=true)
    IEEE519.analyze_csv(csv_filename; fundamental_freq=fundamental_freq, pcc_kv=pcc_kv, plot_results=plot_results)
end

end # module IEEE519API
