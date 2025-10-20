# Start from the official R base image for reproducibility
FROM rocker/r-base:4.3.2

# Install all necessary R packages in one go
# Added tidytext, wordcloud, and its dependency RColorBrewer
RUN R -e 'install.packages(c("data.table", "tidyverse", "tidytext", "wordcloud", "RColorBrewer"), repos="https://cloud.r-project.org")'

# Set a working directory inside the container to keep things tidy
WORKDIR /app

# Copy the R script from your local machine into the container's working directory
# IMPORTANT: The script must be in the same folder as this Dockerfile when you build
COPY BIOS611-HW4.R .

# Set the default command to run when the container starts.
# This will execute your R script.
CMD ["Rscript", "BIOS611-HW4.R"]
