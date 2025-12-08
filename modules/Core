module IEEE519

using FFTW
using Statistics
using Plots
using CSV
using DataFrames

# Define exported functions that users can call
export run_analysis, generate_dummy_csv

# Attempt to initialize plotting backend on module load
function __init__()
    try
        # Try to use PlotlyJS for interactive plots if available
        plotlyjs()
    catch
        # Fallback silently to GR if PlotlyJS isn't fully set up
        gr()
    end
end


# ==========================================
# Data Structures
# ==========================================
struct HarmonicAnalysisResult
    fundamental_freq::Float64
    fundamental_magnitude::Float64
    thd_percent::Float64
    # Using Any for vectors to allow Plotly to handle them easily, though Float64/Int is stricter
    harmonics_magnitudes::Vector{Float64}
    harmonic_orders::Vector{Int}
    freq_vector::Vector{Float64}
    fft_magnitudes::Vector{Float64}
end

# ==========================================
# Core Analysis Functions
# ==========================================

"""
    analyze_harmonics(signal, fs, target_fundamental; max_order=50, enable_window=true)

Performs FFT on time-domain signal to extract harmonic magnitudes and calculate THD.

**Important Note on CSV Data:** Real-world data rarely captures exactly an integer
number of cycles. This causes "spectral leakage" in the FFT. This function applies
a Hanning window by default to mitigate leakage.
"""
function analyze_harmonics(signal::AbstractVector{Float64}, fs::Float64, target_fundamental::Float64; max_order=50, enable_window=true)
    N = length(signal)
    
    processed_signal = signal

    # --- Step A: Apply Windowing ---
    # Crucial for real-world CSV data to prevent spectral leakage
    if enable_window
        # Hanning window calculation
        window = 0.5 .* (1 .- cos.(2π .* (0:(N-1)) ./ (N-1)))
        processed_signal = signal .* window
        # Windowing reduces the energy. We must compensate by a factor.
        # For Hanning, the coherent gain compensation is factor of 2.
        window_compensation = 2.0
    else
        window_compensation = 1.0
    end


    # --- Step B: Perform FFT ---
    F = fft(processed_signal)

    # Calculate magnitudes.
    # 1. Normalize by N for standard FFT definition.
    # 2. Multiply by 2 to get single-sided peak amplitudes for AC components.
    # 3. Multiply by window compensation factor.
    magnitudes_full = (abs.(F) ./ N) .* 2 .* window_compensation

    # We only need the first half of the spectrum (up to Nyquist frequency)
    num_unique_points = floor(Int, N / 2) + 1
    magnitudes = magnitudes_full[1:num_unique_points]
    freqs = range(0, stop=fs/2, length=num_unique_points)


    # --- Step C: Find Fundamental and Harmonics ---
    # Find bin closest to target fundamental (skipping DC at index 1)
    # We look in a range around the target to ensure we don't pick noise
    search_range_idx = (freqs .> target_fundamental * 0.8) .& (freqs .< target_fundamental * 1.2)
    valid_indexes = findall(search_range_idx)
    
    if isempty(valid_indexes)
        error("Could not find a distinct fundamental frequency near $target_fundamental Hz.")
    end

    # Find peak within that search range
    local_peak_idx = argmax(magnitudes[valid_indexes])
    fundamental_idx = valid_indexes[local_peak_idx]

    actual_fundamental_freq = freqs[fundamental_idx]
    fundamental_mag = magnitudes[fundamental_idx]

    harmonic_mags = Float64[]
    harmonic_orders = Int[]
    sum_squares_harmonics = 0.0

    for h in 2:max_order
        target_h_freq = actual_fundamental_freq * h
        # Find closest bin to the target harmonic frequency
        h_idx = argmin(abs.(freqs .- target_h_freq))

        # NOTE ON IEEE 519 COMPLIANCE:
        # This is simple peak-picking. Strict IEEE 519/IEC 61000-4-7 compliance 
        # requires "Harmonic Grouping" (summing energy of adjacent interharmonic bins).
        # This implementation is a simplified approximation.
        if h_idx <= num_unique_points
            h_mag = magnitudes[h_idx]
            push!(harmonic_mags, h_mag)
            push!(harmonic_orders, h)
            sum_squares_harmonics += h_mag^2
        end
    end

    # --- Step D: Calculate THD ---
    # THD = sqrt(sum(h_mag^2)) / fundamental_mag * 100
    thd_pct = fundamental_mag > 1e-9 ? (sqrt(sum_squares_harmonics) / fundamental_mag) * 100.0 : 0.0

    return HarmonicAnalysisResult(
        actual_fundamental_freq, fundamental_mag, thd_pct,
        harmonic_mags, harmonic_orders, freqs, magnitudes
    )
end

function check_voltage_limits(pcc_voltage_kv::Float64, v_thd_pct::Float64)
    limit_vthd = 0.0
    
    # IEEE 519-2014 Table 1 Simplified Limits
    if pcc_voltage_kv <= 1.0
        limit_vthd = 8.0
    elseif pcc_voltage_kv <= 69.0
        limit_vthd = 5.0
    elseif pcc_voltage_kv <= 161.0
        limit_vthd = 2.5
    else
        limit_vthd = 1.5
    end

    status = v_thd_pct <= limit_vthd ? "PASS" : "FAIL"

    println("\n--- IEEE 519 Voltage Compliance Check (Simplified Table 1) ---")
    println("PCC Voltage Level: $pcc_voltage_kv kV")
    println("Measured V_THD: $(round(v_thd_pct, digits=2)) %")
    println("Limit V_THD:    $limit_vthd %")
    println("Status: ** $status **")
end

# ==========================================
# Helper Functions
# ==========================================

"""
    generate_dummy_csv(filename="power_data.csv"; fs=5000.0, duration=0.2, freq=60.0)

Helper function to create synthetic data for testing.
"""
function generate_dummy_csv(filename="power_data.csv"; fs=5000.0, duration=0.2, freq=60.0)
    println("Creating dummy CSV file: $filename")
    t_dummy = 0:(1/fs):duration
    # Create signals with 3rd, 5th, 7th harmonics
    # Voltage: Slightly distorted
    v_dummy = 10000 .* sin.(2π*freq .* t_dummy) .+ 
               500 .* sin.(2π*(3*freq) .* t_dummy) .+ 
               300 .* sin.(2π*(5*freq) .* t_dummy)
    # Current: More heavily distorted (typical of non-linear loads)
    i_dummy = 100 .* sin.(2π*freq .* t_dummy .- 0.5) .+ 
               20 .* sin.(2π*(3*freq) .* t_dummy .- 0.5) .+ 
               15 .* sin.(2π*(5*freq) .* t_dummy .- 0.6) .+
               10 .* sin.(2π*(7*freq) .* t_dummy .- 0.7)

    dummy_df = DataFrame(Time = t_dummy, VoltageA = v_dummy, CurrentA = i_dummy)
    CSV.write(filename, dummy_df)
    println("File created successfully.")
end

# ==========================================
# Main Execution Function
# ==========================================

"""
    run_analysis(csv_filename::String, pcc_kv::Float64, f_sys_target::Float64)

Main entry point. Loads the CSV, calculates Fs, runs harmonic analysis on VoltageA
and CurrentA, checks voltage limits, and plots results.
"""
function run_analysis(csv_filename::String; pcc_kv::Float64=12.47, f_sys_target::Float64=60.0)
    
    if !isfile(csv_filename)
        error("CSV file not found: $csv_filename. Need data to analyze.")
    end

    println("Loading data from: $csv_filename...")
    df = CSV.read(csv_filename, DataFrame)
    
    # Validate Columns
    required_cols = ["Time", "VoltageA", "CurrentA"]
    if !issubset(required_cols, names(df))
        error("CSV missing required columns. Needs: $required_cols. Found: $(names(df))")
    end

    # Extract data so it's clearly typed as Float64 vectors for the analysis function
    t_vec = Float64.(df.Time)
    v_signal = Float64.(df.VoltageA)
    i_signal = Float64.(df.CurrentA)

    # Calculate Sampling Frequency (Fs)
    dt_avg = mean(diff(t_vec))
    Fs_calculated = 1.0 / dt_avg
    println("Data Duration: $(round(t_vec[end] - t_vec[1], digits=3)) s")
    println("Calculated Fs: $(round(Fs_calculated, digits=0)) Hz")

    # Nyquist check
    if Fs_calculated < (50 * f_sys_target * 2)
         @warn "Sampling frequency ($Fs_calculated Hz) might be too low to accurately capture the 50th harmonic (~$(50*f_sys_target) Hz)."
    end

    println("\n--- Starting Harmonic Analysis (Hanning Window Active) ---")

    # Analyze Voltage
    v_results = analyze_harmonics(v_signal, Fs_calculated, f_sys_target; enable_window=true)

    # Analyze Current
    # NOTE: IEEE 519 requires TDD (Total Demand Distortion) for current compliance.
    # TDD uses the maximum historical demand load current ($I_L$) in the denominator.
    # This function calculates THD using the *present* snapshot fundamental current ($I_1$).
    i_results = analyze_harmonics(i_signal, Fs_calculated, f_sys_target; enable_window=true)

    # Display Results
    println("\nVoltage Results (Phase A):")
    println("  Fundamental: $(round(v_results.fundamental_freq, digits=2)) Hz @ $(round(v_results.fundamental_magnitude, digits=1)) Vpk")
    println("  Voltage THD: $(round(v_results.thd_percent, digits=2)) %")

    println("\nCurrent Results (Phase A Snapshot):")
    println("  Fundamental: $(round(i_results.fundamental_freq, digits=2)) Hz @ $(round(i_results.fundamental_magnitude, digits=1)) Apk")
    println("  Current THD: $(round(i_results.thd_percent, digits=2)) % (IMPORTANT: This is snapshot THD, not IEEE 519 TDD)")

    # Check Voltage Limits
    check_voltage_limits(pcc_kv, v_results.thd_percent)

    # --- Visualization ---
    println("\nGenerating Plots...")
    samples_to_plot = min(2000, length(t_vec))

    # 1. Time Domain View
    p1 = plot(t_vec[1:samples_to_plot], v_signal[1:samples_to_plot],
         title="Time Domain Waveforms (Zoomed)", label="Voltage A (V)", ylabel="Amplitude",
         linewidth=1.5)
    # Scale current arbitrarily by 10x just so it visible on the same axis as voltage
    plot!(p1, t_vec[1:samples_to_plot], i_signal[1:samples_to_plot] .* 10,
          label="Current A (x10 A - Scaled)", linewidth=1.5) 

    # 2. Frequency Domain View (Spectrum)
    # Limit X-axis up to 40th harmonic for clarity
    max_freq_plot = f_sys_target * 40

    # Use a log scale for Y-axis to see smaller harmonics better
    p2 = plot(v_results.freq_vector, v_results.fft_magnitudes,
         title="Voltage Frequency Spectrum (FFT + Hanning)",
         xlabel="Frequency (Hz)", ylabel="Magnitude (Vpk)",
         label="Voltage FFT", xlims=(0, max_freq_plot), yscale=:log10)
         
    p3 = bar(v_results.harmonic_orders, v_results.harmonics_magnitudes,
        title="Voltage Harmonic Magnitudes", xlabel="Harmonic Order (h)", ylabel="Vpk",
        label="Vh", xlims=(1.5, 25))

    # Combine plots into a layout
    final_plot = plot(p1, p2, p3, layout=(3,1), size=(900, 1000))
    display(final_plot)
end

end # module IEEE519