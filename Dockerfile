# Use the official Julia image (Debian-based is safer for Plots/GR dependencies)
FROM julia:1.10-bookworm

# 1. Install Linux dependencies required by GR/Plots
# This prevents "library not found" errors on the server
RUN apt-get update && apt-get install -y \
    ca-certificates \
    libgl1-mesa-glx \
    libglib2.0-0 \
    qt6-base-dev \
    xorg \
    && rm -rf /var/lib/apt/lists/*

# 2. Set working directory
WORKDIR /app

# 3. Copy Project files first (for caching layers)
COPY Project.toml Manifest.toml ./

# 4. Install and Precompile Packages
# We run this during BUILD so the startup on Render is faster
# and doesn't crash due to memory limits.
RUN julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'

# 5. Trigger a "warm-up" (Optional but recommended)
# This forces Plots/GR to compile their heavy functions now, not when the user clicks upload.
RUN julia --project=. -e 'using Plots; gr(); plot(rand(5))'

# 6. Copy the rest of the app
COPY . .

# 7. Set Environment Variables
# GKSwstype=100 tells GR to run in "headless" mode (no screen)
ENV GKSwstype=100
ENV JULIA_PROJECT=.
# Restrict threads to save memory on Free Tier
ENV JULIA_NUM_THREADS=1

# 8. Expose the port
EXPOSE 8080

# 9. Start the command
CMD ["julia", "src/IEEE519_API.jl"]