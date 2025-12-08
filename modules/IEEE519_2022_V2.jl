using FFTW
using Statistics
using DataFrames
using Printf
using CSV
using Plots # <--- NEW: For Plotting
gr()  # Use GR backend for plotting

# --- CONFIGURATION STRUCTURE ---
struct IEEE519Config
    f_sys::Float64          # System Frequency (User Input)
    fs::Float64             # Sampling Frequency (Calculated from file)
    I_sc::Float64           # Short Circuit Current
    I_L::Float64            # Max Demand Load Current
    V_nominal::Float64      # Nominal Voltage
end

# --- CORE ANALYSIS FUNCTION ---
function analyze_window(v_win::Vector{Float64}, i_win::Vector{Float64}, cfg::IEEE519Config)
    N = length(v_win)
    
    # FFT
    V_fft = rfft(v_win)
    I_fft = rfft(i_win)
    
    # Normalize to RMS
    mag_V = (abs.(V_fft) ./ N) .* 2 ./ sqrt(2)
    mag_I = (abs.(I_fft) ./ N) .* 2 ./ sqrt(2)
    
    harmonics_V = Float64[]
    harmonics_I = Float64[]
    
    # Harmonic Subgroups (IEC 61000-4-7)
    for h in 1:50
        center_idx = 1 + round(Int, (h * cfg.f_sys) / (cfg.fs / N))
        indices = (center_idx-1):(center_idx+1)
        valid_indices = filter(x -> x > 1 && x <= length(mag_V), indices)
        
        push!(harmonics_V, sqrt(sum(mag_V[valid_indices].^2)))
        push!(harmonics_I, sqrt(sum(mag_I[valid_indices].^2)))
    end
    
    # Calculations
    V_fund = harmonics_V[1]
    I_fund_inst = harmonics_I[1]
    
    V_rss = sqrt(sum(harmonics_V[2:end].^2))
    I_rss = sqrt(sum(harmonics_I[2:end].^2))
    
    THDv = (V_rss / V_fund) * 100
    TDD = (I_rss / cfg.I_L) * 100
    
    return (
        THDv = THDv,
        TDD = TDD,
        I_fund = I_fund_inst,
        Harmonics_V = harmonics_V,
        Harmonics_I = harmonics_I
    )
end

# --- PLOTTING HELPER ---
function generate_plots(v_win, i_win, results, fs, f_sys)
    # Time Axis
    t = (0:length(v_win)-1) ./ fs

    # 1. Time Domain Plots
    p1 = plot(t, v_win, title="Voltage Waveform", ylabel="Volts", label="V", 
              lc=:blue, lw=1.5, xlabel="Time (s)")
    p2 = plot(t, i_win, title="Current Waveform", ylabel="Amps", label="I", 
              lc=:red, lw=1.5, xlabel="Time (s)")

    # 2. Frequency Domain Plots (Bar Charts)
    # Convert to % of Fundamental for visualization
    h_orders = 1:length(results.Harmonics_V)
    
    # We remove DC (index 0 implies strictly 1:50 here)
    # Typically we plot harmonics 2 to 50 for clarity, or 1 to 50.
    # Let's plot 1-25 for better visibility of lower order issues
    orders_view = 1:25
    
    v_spect = (results.Harmonics_V[orders_view] ./ results.Harmonics_V[1]) .* 100
    i_spect = (results.Harmonics_I[orders_view] ./ results.Harmonics_I[1]) .* 100

    p3 = bar(orders_view, v_spect, title="Voltage Spectrum", ylabel="% of Fund", 
             label="V Harmonics", color=:blue, xlabel="Harmonic Order")
    p4 = bar(orders_view, i_spect, title="Current Spectrum", ylabel="% of Fund", 
             label="I Harmonics", color=:red, xlabel="Harmonic Order")

    # Combine into grid
    final_plot = plot(p1, p3, p2, p4, layout=(2,2), size=(1000, 800))
    display(final_plot)
end

# --- INPUT HELPER ---
function get_system_frequency()
    print("\nEnter System Frequency (Hz) [Default: 50.0]: ")
    try
        input_str = strip(readline())
        if isempty(input_str)
            return 50.0
        else
            return parse(Float64, input_str)
        end
    catch
        println("Invalid input. Defaulting to 50.0 Hz")
        return 50.0
    end
end

# --- MAIN RUNNER ---
function run_analysis_on_file(file_path::String)
    println("--- IEEE 519 Analyzer ---")
    println("Loading data from: $file_path")
    
    # 1. Get User Input
    f_sys_user = get_system_frequency()
    println("Using System Frequency: $f_sys_user Hz")

    # 2. Load File
    df = CSV.read(file_path, DataFrame)
    
    # *** MAPPING: ENSURE THESE MATCH YOUR CSV ***
    col_time = "Time"    
    col_volt = "Voltage" 
    col_curr = "Current" 
    
    if !all(c -> c in names(df), [col_time, col_volt, col_curr])
        println("ERROR: Column names mismatch. Found: ", names(df))
        return
    end

    t_raw = df[!, col_time]
    volts = df[!, col_volt]
    amps = df[!, col_curr]
    
    # 3. Calc Sampling Freq
    dt = mean(diff(t_raw))
    fs_calculated = 1.0 / dt
    @printf("Calculated Sampling Freq: %.1f Hz\n", fs_calculated)

    config = IEEE519Config(
        f_sys_user,      # <--- User Input
        fs_calculated,   
        2000.0,          # I_sc
        100.0,           # I_L
        480.0            # V_nominal
    )
    
    # 4. Select Window (Approx 12 cycles)
    cycles_to_analyze = 12
    samples_needed = round(Int, (cycles_to_analyze / f_sys_user) * fs_calculated)
    
    if length(volts) < samples_needed
        println("Warning: Short file. Using available data.")
        samples_needed = length(volts)
    end
    
    v_window = volts[1:samples_needed]
    i_window = amps[1:samples_needed]
    
    # 5. Analyze & Plot
    results = analyze_window(v_window, i_window, config)
    
    # Print Text Report
    println("\n--- Results ---")
    @printf("Voltage THD: %.2f%%\n", results.THDv)
    @printf("Current TDD: %.2f%%\n", results.TDD)
    
    # Generate Plots
    println("Generating plots...")
    generate_plots(v_window, i_window, results, fs_calculated, f_sys_user)
end