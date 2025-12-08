using FFTW
using Statistics
using DataFrames
using Printf
using Pkg

# --- CONFIGURATION STRUCTURE ---
struct IEEE519Config
    f_sys::Float64          # System Frequency (60.0 or 50.0 Hz)
    fs::Float64             # Sampling Frequency (Hz)
    I_sc::Float64           # Short Circuit Current at PCC (Amps)
    I_L::Float64            # Max Demand Load Current (Amps) - 12-month rolling avg peak
    V_nominal::Float64      # Nominal Voltage (Volts)
end

# --- CORE ANALYSIS FUNCTION ---
"""
    analyze_window(voltage_window, current_window, config)

Performs IEC 61000-4-7 Harmonic Subgroup analysis on a single 10/12 cycle window.
Returns a NamedTuple with THDv, TDD, and individual harmonics.
"""
function analyze_window(v_win::Vector{Float64}, i_win::Vector{Float64}, cfg::IEEE519Config)
    # 1. Determine Window Parameters
    # IEEE 519 requires approx 200ms window (10 cycles @ 50Hz, 12 cycles @ 60Hz)
    # This results in exactly 5Hz bin spacing in the FFT.
    N = length(v_win)
    
    # 2. Perform FFT
    # rfft computes real-input FFT (output is N/2 + 1 complex values)
    V_fft = rfft(v_win)
    I_fft = rfft(i_win)
    
    # Normalize FFT magnitude to RMS amplitude
    # (Divide by N, multiply by sqrt(2) for RMS? 
    #  Standard approach: Magnitude = abs(FFT)/N * 2 (for peak). RMS = Peak / sqrt(2).)
    freqs = rfftfreq(N, cfg.fs)
    mag_V = (abs.(V_fft) ./ N) .* 2 ./ sqrt(2)
    mag_I = (abs.(I_fft) ./ N) .* 2 ./ sqrt(2)
    
    # 3. Calculate Harmonic Subgroups (IEC 61000-4-7)
    # Subgroup H_n includes the harmonic bin (C_k) plus adjacent bins (C_k-1, C_k+1)
    # Bin resolution is approx 5Hz.
    # Fundamental bin index (approximate):
    bin_fund = findmin(abs.(freqs .- cfg.f_sys))[2]
    
    harmonics_V = Float64[]
    harmonics_I = Float64[]
    
    # Loop through harmonics up to 50th
    for h in 1:50
        # Calculate expected index for harmonic h
        center_idx = 1 + round(Int, (h * cfg.f_sys) / (cfg.fs / N))
        
        # IEC Grouping: sqrt(C_k^2 + C_k-1^2 + C_k+1^2)
        # We need boundary checks for indices
        indices = (center_idx-1):(center_idx+1)
        valid_indices = filter(x -> x > 1 && x <= length(mag_V), indices)
        
        # Voltage Subgroup
        v_sg = sqrt(sum(mag_V[valid_indices].^2))
        push!(harmonics_V, v_sg)
        
        # Current Subgroup
        i_sg = sqrt(sum(mag_I[valid_indices].^2))
        push!(harmonics_I, i_sg)
    end
    
    # 4. Calculate THD (Voltage) and TDD (Current)
    # THDv uses V_fundamental (h=1) as denominator
    # TDD uses I_L (Max Demand Load Current) as denominator
    
    V_fund = harmonics_V[1]
    I_fund_inst = harmonics_I[1] # Instantaneous fundamental (not used for TDD denominator)
    
    # RSS of harmonics 2 to 50
    V_rss = sqrt(sum(harmonics_V[2:end].^2))
    I_rss = sqrt(sum(harmonics_I[2:end].^2))
    
    THDv = (V_rss /cfg.V_nominal) * 100  # Usually % of Nominal or Fundamental. IEEE uses Nominal for limits? 
                                      # Actually IEEE defines THD as % of Fundamental for Voltage.
    THDv_calc = (V_rss / V_fund) * 100
    
    TDD = (I_rss / cfg.I_L) * 100     # CURRENT uses I_L (Max Demand)
    
    return (
        THDv = THDv_calc,
        TDD = TDD,
        I_fund = I_fund_inst,
        Harmonics_V = harmonics_V,
        Harmonics_I = harmonics_I
    )
end

# --- LIMIT CHECKING FUNCTION ---
function check_compliance(results, cfg::IEEE519Config)
    SCR = cfg.I_sc / cfg.I_L
    println("\n--- IEEE 519-2022 Compliance Report ---")
    @printf("System: %.1f Hz | SCR: %.1f\n", cfg.f_sys, SCR)
    @printf("Measured V_fund: %.2f V | Measured I_fund: %.2f A\n", results.Harmonics_V[1], results.I_fund)
    println("-----------------------------------------")

    # 1. Voltage Check (Table 1)
    # Limits depend on Voltage Level (assuming <= 69kV for this example)
    # Limit: 8.0% for <= 1kV, 5.0% for 1-69kV
    v_limit_thd = cfg.V_nominal <= 1000 ? 8.0 : 5.0
    v_status = results.THDv <= v_limit_thd ? "PASS" : "FAIL"
    
    @printf("Voltage THD: %.2f%% (Limit: %.1f%%) -> %s\n", results.THDv, v_limit_thd, v_status)

    # 2. Current Check (Table 2 - simplified for < 69kV)
    # Determine TDD limit based on SCR
    if SCR < 20
        tdd_limit = 5.0
    elseif SCR < 50
        tdd_limit = 8.0
    elseif SCR < 100
        tdd_limit = 12.0
    elseif SCR < 1000
        tdd_limit = 15.0
    else
        tdd_limit = 20.0
    end
    
    i_status = results.TDD <= tdd_limit ? "PASS" : "FAIL"
    @printf("Current TDD: %.2f%% (Limit: %.1f%% @ SCR %.0f) -> %s\n", results.TDD, tdd_limit, SCR, i_status)
    
    # 3. Individual Harmonic Check (Example for 5th Harmonic)
    # IEEE 519 2022 allows Even harmonics relaxed limits.
    # This is a basic check for TDD-level compliance.
    
    return (THDv_Pass = (results.THDv <= v_limit_thd), TDD_Pass = (results.TDD <= tdd_limit))
end

# --- DEMONSTRATION ---
function run_demo()
    # Mock Data Generation
    # Create a 60Hz system with sampling at 3072Hz (51.2 samples/cycle)
    fs = 3072.0
    f_sys = 60.0
    samples_per_window = round(Int, fs / f_sys * 12) # 12 cycle window for 60Hz
    t = (0:samples_per_window-1) / fs
    
    # Generate Voltage: Clean 480V with slight 5th harmonic
    v_nominal = 480.0
    volts = (v_nominal * sqrt(2)) .* sin.(2π * f_sys * t) .+ 
            (v_nominal * 0.04 * sqrt(2)) .* sin.(2π * 5 * f_sys * t) # 4% 5th harmonic
            
    # Generate Current: Distorted Load (e.g., VFD)
    # High 5th (20%) and 7th (10%)
    i_max_demand = 100.0 # Our facility is rated for 100A
    i_fund = 80.0 # Currently running at 80A
    amps = (i_fund * sqrt(2)) .* sin.(2π * f_sys * t) .+ 
           (i_fund * 0.20 * sqrt(2)) .* sin.(2π * 5 * f_sys * t) .+
           (i_fund * 0.10 * sqrt(2)) .* sin.(2π * 7 * f_sys * t)
           
    # Configuration
    config = IEEE519Config(
        60.0,   # f_sys
        fs,     # fs
        2000.0, # I_sc (Short circuit current - utility provided)
        i_max_demand,  # I_L (Max demand load current)
        v_nominal # V_nominal
    )
    
    # Analyze
    results = analyze_window(volts, amps, config)
    check_compliance(results, config)
end

# Run the demo
run_demo() 