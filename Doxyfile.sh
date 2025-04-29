#-------------------------------------------------------------
# Doxygen basic configuration file for C++ project
#-------------------------------------------------------------

# Input directory containing your .cpp and .h files
INPUT = ./src

# Output directory for documentation
OUTPUT_DIRECTORY = ./docs

# File patterns to consider (C++ files)
FILE_PATTERNS = *.cpp *.h

# Enable recursive search into subfolders
RECURSIVE = YES

# Extract documentation from comments
EXTRACT_ALL = YES

# Generate HTML documentation
GENERATE_HTML = YES

# (Optional) Also generate LaTeX documentation
GENERATE_LATEX = NO

# Use Markdown support (if you want to write readme.md files)
MARKDOWN_SUPPORT = YES

# Show undocumented members in output (for full completeness)
SHOW_UNDOCUMENTED = YES

# Show warnings for undocumented members
WARN_IF_UNDOCUMENTED = YES

# Treat struct comments and class comments equally
EXTRACT_STRUCTS = YES
EXTRACT_CLASSES = YES

# Sort members by appearance in source file (easier to read)
SORT_MEMBER_DOCS = NO

# Strip path from file names in documentation
STRIP_FROM_PATH = ./src