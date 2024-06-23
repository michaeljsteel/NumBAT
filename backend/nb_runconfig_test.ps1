# Define output file path
$OUTPUT_FILE = Join-Path (Get-Location) -ChildPath "nb_build_config.txt"

# Check if backend directory exists (assuming NumBAT root)
if (!(Test-Path -Path "backend")) {
  Write-Host "This script must be executed from the NumBAT root directory."
  exit 1
}

# Clear the output file
Clear-Content -Path $OUTPUT_FILE

# Write header information
Add-Content -Path $OUTPUT_FILE -Value "NumBAT config test file"
Add-Content -Path $OUTPUT_FILE -Value "------------------------"

# Get system information
Add-Content -Path $OUTPUT_FILE -Value (Get-Date -Format "Date:        {0}" )

# Limited OS information retrieval on Windows
Add-Content -Path $OUTPUT_FILE -Value (Get-WmiObject Win32_OperatingSystem | Format-List @{Name="Env";Expression={$_.caption + " " + $_.version}})

# LD_LIBRARY_PATH not used directly on Windows
Write-Host "Windows doesn't use LD_LIBRARY_PATH."

# Find Python and Python3 paths (might need adjustment based on installation)
$python_path = Where-Object -Path python -ListDirectory

# Check for Python versions (if Python paths are valid)
if ($python_path) {
  $python_version = & python -Command "& { import sys; print(sys.version) }"
  Add-Content -Path $OUTPUT_FILE -Value ("python path:  " + $python_path)
  Add-Content -Path $OUTPUT_FILE -Value ("python ver.:  " + $python_version)
}

$python3_path = Where-Object -Path python3 -ListDirectory

# Check for Python3 versions (if Python paths are valid)
if ($python3_path) {
  $python3_version = & python3 -Command "& { import sys; print(sys.version) }"
  Add-Content -Path $OUTPUT_FILE -Value ("python3 path: " + $python3_path)
  Add-Content -Path $OUTPUT_FILE -Value ("python3 ver.: " + $python3_version)
}

# Check for numpy, scipy, matplotlib (if Python paths are valid)
$libraries = @("numpy", "scipy", "matplotlib")
foreach ($library in $libraries) {
  $version_cmd = "& { import $library; print($library.__version__) }"
  try {
    $version = if ($python_path) { & python -Command $version_cmd } else { & python3 -Command $version_cmd }
    Add-Content -Path $OUTPUT_FILE -Value ("$library ver.:  " + $version)
  } catch {
    Add-Content -Path $OUTPUT_FILE -Value ("$library not found")
  }
}

# Write NumBAT build settings
Add-Content -Path $OUTPUT_FILE -Value ""
Add-Content -Path $OUTPUT_FILE -Value "NumBAT build settings:"
Add-Content -Path $OUTPUT_FILE -Value ""
Get-Content -Path "backend/fortran/nb_build_config.txt" | Add-Content -Path $OUTPUT_FILE

Write-Host ("Done! Information written to '{0}'.") -f $OUTPUT_FILE

