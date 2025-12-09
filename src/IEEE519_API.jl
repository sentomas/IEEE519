using Oxygen
using HTTP
using FFTW
using Statistics
using DataFrames
using CSV
using Plots
using Base64
using JSON3

# --- 1. CONFIGURATION & SETUP ---
ENV["GKSwstype"] = "100" 
gr()

struct IEEE519Config
    f_sys::Float64
    fs::Float64
    I_sc::Float64
    I_L::Float64
    V_nominal::Float64
end

# --- 2. CORE LOGIC ---
function analyze_window(v_win::Vector{Float64}, i_win::Vector{Float64}, cfg::IEEE519Config)
    N = length(v_win)
    
    V_fft = rfft(v_win)
    I_fft = rfft(i_win)
    
    mag_V = (abs.(V_fft) ./ N) .* 2 ./ sqrt(2)
    mag_I = (abs.(I_fft) ./ N) .* 2 ./ sqrt(2)
    
    harmonics_V = Float64[]
    harmonics_I = Float64[]
    
    for h in 1:50
        center_idx = 1 + round(Int, (h * cfg.f_sys) / (cfg.fs / N))
        indices = (center_idx-1):(center_idx+1)
        valid_indices = filter(x -> x > 1 && x <= length(mag_V), indices)
        
        push!(harmonics_V, sqrt(sum(mag_V[valid_indices].^2)))
        push!(harmonics_I, sqrt(sum(mag_I[valid_indices].^2)))
    end
    
    V_fund = harmonics_V[1]
    I_fund_inst = harmonics_I[1]
    
    V_rss = sqrt(sum(harmonics_V[2:end].^2))
    I_rss = sqrt(sum(harmonics_I[2:end].^2))
    
    THDv = (V_rss / V_fund) * 100
    TDD = (I_rss / cfg.I_L) * 100
    
    return (
        THDv = THDv,
        TDD = TDD,
        Harmonics_V = harmonics_V,
        Harmonics_I = harmonics_I
    )
end

# --- 3. PLOTTING ---

function generate_plot_base64(v_win, i_win, results, fs)
    t = (0:length(v_win)-1) ./ fs
    
    p1 = plot(t, v_win, title="Voltage", ylabel="V", label="", lc=:blue, titlefontsize=10)
    p2 = plot(t, i_win, title="Current", ylabel="A", label="", lc=:red, titlefontsize=10)
    
    orders = 1:25
    v_spect = (results.Harmonics_V[orders] ./ results.Harmonics_V[1]) .* 100
    i_spect = (results.Harmonics_I[orders] ./ results.Harmonics_I[1]) .* 100
    
    p3 = bar(orders, v_spect, title="V Harmonics %", label="", color=:blue, titlefontsize=10)
    p4 = bar(orders, i_spect, title="I Harmonics %", label="", color=:red, titlefontsize=10)
    
    # Corrected margin syntax (5*Plots.mm)
    final_plot = plot(p1, p3, p2, p4, layout=(2,2), size=(800, 600), margin=5*Plots.mm)
    
    io = IOBuffer()
    
    # --- FIX IS HERE ---
    # We use 'show' with a MIME type instead of 'savefig'
    show(io, MIME("image/png"), final_plot) 
    # -------------------
    
    return base64encode(take!(io))
end

# --- 4. THE API ENDPOINT ---

@post "/analyze" function(req::HTTP.Request)
    try
        # A. Parse Form Data
        form_data = HTTP.parse_multipart_form(req)
        
        file_entry = nothing
        f_sys_val = 50.0 
        
        for part in form_data
            if part.name == "csv_file"
                file_entry = part
            elseif part.name == "frequency"
                f_sys_val = parse(Float64, String(read(part.data)))
            end
        end

        if file_entry === nothing
            return HTTP.Response(400, "Error: No file uploaded.")
        end

        # --- B. ROBUST CSV READING (THE FIX) ---
        # We use read() to extract raw bytes, ignoring whether it's a Stream or Buffer
        raw_bytes = read(file_entry.data)
        df = CSV.read(IOBuffer(raw_bytes), DataFrame)
        
        # --- DEBUG: Print columns to help you if it fails ---
        println("Processed CSV Columns: ", names(df))

        # --- C. MAPPING (ENSURE THESE MATCH YOUR FILE) ---
        # If your file has "Time_s", change "Time" to "Time_s" below
        t_raw = df[!, "Time"]
        volts = df[!, "Voltage"]
        amps = df[!, "Current"]
        
        # D. Calculate FS
        dt = mean(diff(t_raw))
        fs_calculated = 1.0 / dt
        
        # E. Config
        cfg = IEEE519Config(f_sys_val, fs_calculated, 2000.0, 100.0, 480.0)
        
        # F. Select Window
        samples_needed = round(Int, (12 / f_sys_val) * fs_calculated)
        len = min(length(volts), samples_needed)
        v_win = volts[1:len]
        i_win = amps[1:len]
        
        # G. Analyze & Plot
        results = analyze_window(v_win, i_win, cfg)
        img_base64 = generate_plot_base64(v_win, i_win, results, fs_calculated)
        
        # H. Response
        response_data = Dict(
            "status" => "success",
            "system_freq" => f_sys_val,
            "sampling_freq" => round(fs_calculated, digits=1),
            "thd_v" => round(results.THDv, digits=2),
            "tdd_i" => round(results.TDD, digits=2),
            "harmonics_v" => results.Harmonics_V[1:25],
            "plot_image" => "data:image/png;base64,$img_base64"
        )
        
        return json(response_data)
        
    catch e
        println("\n!!! SERVER ERROR !!!")
        showerror(stdout, e, catch_backtrace())
        return HTTP.Response(500, "Server Error: $(e)")
    end
end
# --- ADD THIS TO JULIA FILE FOR HEALTH CHECK ---
# 1. Serve the Dashboard at the root URL
@get "/" function()
    return file("src/index.html")
end

# 2. Move the Health Check to a specific path (Optional but good)
@get "/health" function()
    return "API Online"
end
# --- 5. SERVER SETUP WITH CORS ---

function CorsMiddleware(handler)
    return function(req::HTTP.Request)
        if req.method == "OPTIONS"
            return HTTP.Response(200, [
                "Access-Control-Allow-Origin" => "*",
                "Access-Control-Allow-Methods" => "POST, GET, OPTIONS",
                "Access-Control-Allow-Headers" => "*"
            ])
        end
        res = handler(req)
        HTTP.setheader(res, "Access-Control-Allow-Origin" => "*")
        return res
    end
end

println("Starting IEEE 519 API on 0.0.0.0:8080")
serve(host="0.0.0.0", port=8080, middleware=[CorsMiddleware])