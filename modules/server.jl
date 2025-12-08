using HTTP, JSON3
using CSV, DataFrames
include("IEEE519.jl")
include("api.jl")

function parse_multipart(body::String, boundary::String)
    """Parse multipart form data"""
    parts = Dict{String, Any}()
    sections = split(body, "--$boundary")
    for section in sections
        if isempty(strip(section)) || section == "--"
            continue
        end
        lines = split(section, "\r\n")
        headers_done = false
        current_key = ""
        content = []
        for line in lines
            if !headers_done
                if contains(line, "name=")
                    # Extract name from Content-Disposition header
                    m = match(r"name=\"([^\"]+)\"", line)
                    if m !== nothing
                        current_key = m.captures[1]
                    end
                end
                if line == ""
                    headers_done = true
                end
            else
                push!(content, line)
            end
        end
        if !isempty(current_key)
            parts[current_key] = join(content[1:end-1], "\r\n")
        end
    end
    parts
end

function handle_request(req::HTTP.Request)
    # ... (OPTIONS check) ...
    
    # This is the ONLY route defined for processing data
    if req.method == "POST" && req.target == "/analyze"
        # ... (logic to process uploaded file) ...
    end
    
    # For any other request (like GET / from a browser), return 404 Not Found
    return HTTP.Response(404, "Not found")
end
    
    if req.method == "POST" && req.target == "/analyze"
        try
            ct = HTTP.header(req, "Content-Type", "")
            body = String(req.body)
            
            # Parse boundary from Content-Type
            boundary = match(r"boundary=([^\s;]+)", ct)
            if boundary === nothing
                return HTTP.Response(400, "No boundary in Content-Type")
            end
            boundary = boundary.captures[1]
            
            parts = parse_multipart(body, boundary)
            if !haskey(parts, "file") || !haskey(parts, "fundamental_freq") || !haskey(parts, "pcc_kv")
                return HTTP.Response(400, "Missing required form fields")
            end
            
            # Save uploaded CSV
            filename = "uploaded_$(Int(time())).csv"
            open(filename, "w") do f
                write(f, parts["file"])
            end
            
            freq = parse(Float64, parts["fundamental_freq"])
            pcc = parse(Float64, parts["pcc_kv"])
            
            # Call analysis
            result = IEEE519API.analyze_csv(filename; fundamental_freq=freq, pcc_kv=pcc, plot_results=false)
            
            # Convert result to JSON-serializable dict
            result_dict = Dict(
                "voltage" => Dict(
                    "fundamental_freq" => result.voltage.fundamental_freq,
                    "fundamental_magnitude" => result.voltage.fundamental_magnitude,
                    "thd_percent" => result.voltage.thd_percent,
                    "harmonics_magnitudes" => result.voltage.harmonics_magnitudes,
                    "harmonic_orders" => result.voltage.harmonic_orders
                ),
                "current" => Dict(
                    "fundamental_freq" => result.current.fundamental_freq,
                    "fundamental_magnitude" => result.current.fundamental_magnitude,
                    "thd_percent" => result.current.thd_percent,
                    "harmonics_magnitudes" => result.current.harmonics_magnitudes,
                    "harmonic_orders" => result.current.harmonic_orders
                ),
                "sampling_freq" => result.sampling_freq,
                "nyquist_ok" => result.nyquist_ok,
                "compliance" => Dict(
                    "limit" => result.compliance.limit,
                    "measured" => result.compliance.measured,
                    "status" => result.compliance.status
                )
            )
            
            return HTTP.Response(200,
                ["Access-Control-Allow-Origin" => "*",
                 "Content-Type" => "application/json"],
                body=JSON3.write(result_dict))
        catch e
            println("Error: $e")
            return HTTP.Response(500,
                ["Access-Control-Allow-Origin" => "*"],
                body="Error: $(string(e))")
        end
    end
    
    return HTTP.Response(404, "Not found")
end

println("Starting IEEE519 API server on http://localhost:8000")
HTTP.serve(handle_request, HTTP.Sockets.localhost, 8000)
