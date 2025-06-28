# Detect OS
OS_NAME := $(shell uname -s)

# Define OS-specific commands
ifeq ($(OS_NAME),Linux)
    CP = cp
    RM = rm -f
    MKDIR = mkdir -p
    # ... other Linux commands
else ifeq ($(OS_NAME),Darwin) # macOS
    CP = cp
    RM = rm -f
    MKDIR = mkdir -p
    # ... other macOS commands
else ifeq ($(OS_NAME),MINGW32_NT) # Git Bash / MinGW
    CP = cp
    RM = rm -f
    MKDIR = mkdir -p
    # ... other Git Bash / MinGW commands
else ifeq ($(OS_NAME),MSYS_NT) # MSYS2
    CP = cp
    RM = rm -f
    MKDIR = mkdir -p
    # ... other MSYS2 commands
else # Assume Windows cmd.exe or PowerShell if uname is not available or recognized
    CP = copy /Y
    RM = del /Q /F
    MKDIR = mkdir
    # You might need to add specific logic for PowerShell if it's your primary target
    # e.g., CP = Copy-Item -Force -Path
endif