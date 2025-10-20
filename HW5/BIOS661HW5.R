# -----------------------------------------------------------------------------
# Computational Homework: Task 1 - Clusters and the Gap Statistic
#
# Description: This script implements a function to generate clustered data.
# It then runs a simulation to test how the Gap Statistic's ability to identify
# the correct number of clusters changes as the clusters move closer together
# across various dimensions.
# -----------------------------------------------------------------------------


library("cluster")
library("plotly")
library("ggplot2")
library("dplyr")
library("tibble")
library("RSpectra")
library("igraph")

# Set a seed for reproducibility.
set.seed(10142025)


# -- 2. DATA GENERATION FUNCTION --
# ---------------------------------
# This function generates 'n' clusters of 'k' points each in an 'n'-dimensional
# space. The cluster centers are located at the positive corners of a hypercube.

generate_hypercube_clusters <- function(n, k, side_length, noise_sd = 1.0) {
  # Validate input
  if (n < 2) stop("Number of dimensions 'n' must be 2 or greater.")
  
  # This will store the points for all clusters
  all_points <- list()
  
  # Loop to create each of the 'n' clusters
  for (i in 1:n) {
    # Create the center for the current cluster.
    # It's a vector of zeros, with the 'side_length' (L) at the i-th position.
    center <- rep(0, n)
    center[i] <- side_length
    
    # Generate 'k * n' random numbers from a Gaussian (normal) distribution.
    # This represents the random scatter of points around the cluster center.
    noise <- matrix(rnorm(k * n, mean = 0, sd = noise_sd),
                    nrow = k,
                    ncol = n)
    
    # Add the center coordinates to the noise matrix.
    cluster_points <- sweep(noise, 2, center, FUN = "+")
    
    # Store the generated points
    all_points[[i]] <- cluster_points
  }
  
  # Combine the lists of points into a single matrix and return it.
  # We also add true cluster labels for potential future validation, though
  # the gap statistic itself does not use them.
  final_data <- do.call(rbind, all_points)
  true_labels <- rep(1:n, each = k)
  
  return(list(data = final_data, labels = true_labels))
}


# -- 3. SIMULATION SETUP & EXECUTION --
# -------------------------------------
# Define the parameters for the simulation as specified in the assignment.
dims_to_test <- c(6, 5, 4, 3, 2)
side_lengths_to_test <- seq(10, 1, by = -1)
K_POINTS <- 100       # Number of points per cluster (k)
NOISE_SD <- 1.0      # Standard deviation of cluster spread

# Create an empty data frame to store the results of our simulation.
simulation_results <- tibble(
  n_dims = integer(),
  side_length = numeric(),
  estimated_k = integer()
)


# --- RUN THE SIMULATION ---
# This set of nested loops will take a few minutes to run.
message("Starting simulation. This may take several minutes...")

for (n in dims_to_test) {
  for (L in side_lengths_to_test) {
    
    # Announce progress to the console
    message(paste0("Running: ", n, " dimensions, side length = ", L))
    
    # 1. Generate the data for the current set of parameters
    generated_data <- generate_hypercube_clusters(n, K_POINTS, L, NOISE_SD)$data
    
    # 2. Run the Gap Statistic analysis
    # K.max is set to n+5 to allow the algorithm to check for more clusters
    # than the true number, which is necessary for it to work correctly.
    # We pass nstart and iter.max directly to clusGap, which will forward
    # them to the kmeans function.
    gap_stat <- clusGap(
      generated_data,
      FUNcluster = kmeans,
      K.max = 10,
      # B = 50, # Number of bootstrap samples.
      nstart = 20,
      iter.max = 50,
      verbose = FALSE # Suppress clusGap progress bar for a cleaner log
    )
    
    # 3. Determine the optimal number of clusters (k)
    # The maxSE method with 'firstSEmax' is a standard, robust way to select k
    # from the gap statistic results. It finds the smallest k for which the
    # gap value is within one standard error of the next k's gap value.
    estimated_k <- maxSE(f = gap_stat$Tab[, "gap"], 
                         SE.f = gap_stat$Tab[, "SE.sim"], 
                         method = "firstSEmax")
    
    # 4. Store the results
    simulation_results <- simulation_results %>%
      add_row(n_dims = n, side_length = L, estimated_k = estimated_k)
  }
}

message("Simulation complete.")


# -- 4. VISUALIZATION OF RESULTS --
# ---------------------------------

# Arrange the results for plotting
plot_data <- simulation_results %>%
  mutate(n_dims_factor = factor(n_dims,
                                levels = dims_to_test,
                                labels = paste(dims_to_test, "Dimensions")))

final_plot <- ggplot(plot_data, aes(x = side_length, y = estimated_k)) +
  geom_hline(aes(yintercept = n_dims), color = "red", linetype = "dashed", linewidth = 1) +
  geom_line(color = "dodgerblue", linewidth = 1) +
  geom_point(color = "dodgerblue", size = 3) +
  facet_wrap(~ n_dims_factor, ncol = 3) +
  scale_x_reverse(breaks = seq(10, 1, by = -1)) +
  scale_y_continuous(breaks = function(limits) seq(floor(limits[1]), ceiling(limits[2]), by = 1)) +
  labs(
    title = "Gap Statistic's Estimate of Cluster Count",
    subtitle = "The red dashed line indicates the true number of clusters.",
    x = "Side Length (L)",
    y = "Estimated Number of Clusters (k)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    strip.text = element_text(face = "bold")
  )

# Print the plot to the RStudio plot pane
print(final_plot)

# ggsave("C:/Eneda_Research/FOXA1_ILC_RNA_DMSOE2/Kmeans_performance.pdf",
#        plot = final_plot,
#        width = 12, height = 8, dpi = 300)
# 
# message("Plot saved")


# -- 5. INTERPRETATION OF RESULTS --
# ----------------------------------
# The plots display the estimated number of clusters (blue line) versus the
# side length (L), which controls cluster separation. The red dashed line
# indicates the true number of clusters (which is equal to 'n', the number
# of dimensions in each facet). The x-axis is reversed to show clusters
# moving from well-separated (L=10) to highly overlapping (L=1).
#
# General Observation:
# For every dimension tested, the Gap Statistic correctly identifies the true
# number of clusters when the side length L is large, indicating clear
# separation. As L decreases, a distinct "breaking point" is reached where
# the clusters become too close for the algorithm to distinguish, and the
# estimated number of clusters abruptly drops to 1.
#
# The Effect of Dimensionality:
# The results show that the clustering method performs WORSE in higher
# dimensions. The "breaking point"—the value of L at which the estimate
# drops—occurs at larger side lengths for higher dimensions. This phenomenon 
# is an example of the "Curse of Dimensionality," where the increasing 
# volume of the space makes density-based estimation (which is implicitly 
# what k-means relies on) more difficult.





# -----------------------------------------------------------------------------
# Computational Homework: Task 2 - Spectral Clustering on Concentric Shells
#
# Description: This script generates 3D concentric shell data, implements a
# spectral clustering algorithm, and uses the Gap Statistic to determine
# at what point the shells become too close for the algorithm to distinguish.
# -----------------------------------------------------------------------------


# -- 1. DATA GENERATION FUNCTION --
# ---------------------------------
# This function generates 'n_shells' concentric shells in 3D space.

# Set a seed for reproducibility.
set.seed(10142025)

generate_shell_clusters <- function(n_shells, k_per_shell, max_radius, noise_sd = 0.1) {
  # Define the mean radius for each shell, evenly spaced from a small
  # inner radius up to the max_radius.
  shell_radii <- seq(from = max_radius / n_shells, to = max_radius, length.out = n_shells)
  
  all_points_list <- lapply(shell_radii, function(r) {
    # Add Gaussian noise to the radius to give the shell thickness
    radii_with_noise <- r + rnorm(k_per_shell, mean = 0, sd = noise_sd)
    
    # Generate points uniformly on a sphere using spherical coordinates
    phi <- runif(k_per_shell, 0, 2 * pi)      # Azimuthal angle
    theta <- acos(runif(k_per_shell, -1, 1)) # Polar angle
    
    # Convert spherical to Cartesian coordinates
    x <- radii_with_noise * sin(theta) * cos(phi)
    y <- radii_with_noise * sin(theta) * sin(phi)
    z <- radii_with_noise * cos(theta)
    
    tibble(x = x, y = y, z = z)
  })
  
  # Add true cluster labels and combine into a single data frame
  all_points <- bind_rows(all_points_list, .id = "true_label") %>%
    mutate(true_label = as.factor(true_label))
  
  return(all_points)
}


# -- 2. VISUALIZE THE DATA STRUCTURE (Interactive Plot) --
# -------------------------------------------------------
# Before running the full simulation, let's generate one sample dataset
# and visualize it to confirm our generator function works as expected.

message("Generating sample data for 3D visualization...")
sample_data <- generate_shell_clusters(n_shells = 4, k_per_shell = 100, max_radius = 10)

# Create an interactive 3D scatter plot
interactive_plot <- plot_ly(
  data = sample_data,
  x = ~x, y = ~y, z = ~z,
  color = ~true_label,
  type = "scatter3d",
  mode = "markers",
  marker = list(size = 3)
) %>% layout(
  title = "Interactive 3D Visualization of Concentric Shells",
  scene = list(xaxis = list(title = 'X'),
               yaxis = list(title = 'Y'),
               zaxis = list(title = 'Z'))
)

# Show the plot. In RStudio, this will appear in the Viewer pane.
print(interactive_plot)
message("Interactive plot created. Check the Viewer pane.")
#htmlwidgets::saveWidget(interactive_plot, "my_3d_plot.html")


# -- 3. SPECTRAL CLUSTERING WRAPPER FUNCTION --
# ---------------------------------------------
# This function wraps the entire spectral clustering logic so it can be
# passed directly to `clusGap`.

# 需要：install.packages(c("igraph", "RSpectra"))
spectral_wrapper <- function(x, k, d_threshold = 1) {
  X <- as.matrix(x)
  n <- nrow(X)
  if (k <= 1L || n <= 1L) return(list(cluster = rep(1L, n)))
  
  # 1) Build a Similarity Graph (Adjacency Matrix A )
  Dmat <- as.matrix(dist(X, method = "euclidean"))
  A <- (Dmat <= d_threshold) * 1L
  diag(A) <- 0L
  A <- (A | t(A)) * 1L
  
  if (sum(A) == 0L) {
    km_fb <- kmeans(scale(X), centers = k, nstart = 10)
    return(list(cluster = as.integer(km_fb$cluster)))
  }
  
  # 2) Compute the Laplacian Matrix ( L )
  g <- igraph::graph_from_adjacency_matrix(A, mode = "undirected", diag = FALSE)
  Lsym <- igraph::laplacian_matrix(g, normalized = TRUE, sparse = TRUE)
  
  # 3) Eigen-decomposition
  kk_req <- min(as.integer(k), n - 1L)

  eig <- RSpectra::eigs_sym(Lsym, k = kk_req, which = "SM")
  U <- eig$vectors
  if (is.null(dim(U))) {
    U <- matrix(U, nrow = n, ncol = 1L) 
  } else {
    U <- as.matrix(U)
  }
  kk_eff <- ncol(U)
  
  if (kk_eff == 0L) {
    km_fb <- kmeans(scale(X), centers = k, nstart = 25)
    return(list(cluster = as.integer(km_fb$cluster)))
  }
  
  # 4) Cluster the Eigenvectors
  rn <- sqrt(rowSums(U^2))
  rn[rn == 0] <- 1
  U_norm <- U / rn

  uniq_count <- nrow(unique(round(U_norm, 8)))
  if (uniq_count < k) {
    jitter_sd <- 1e-6
    U_norm <- U_norm + matrix(rnorm(n * kk_eff, sd = jitter_sd), nrow = n, ncol = kk_eff)
    uniq_count <- nrow(unique(round(U_norm, 8)))
  }
  if (uniq_count < k) {
    km_fb <- kmeans(scale(X), centers = k, nstart = 25)
    return(list(cluster = as.integer(km_fb$cluster)))
  }
  
  res <- tryCatch(
    {
      km <- kmeans(U_norm, centers = k, nstart = 25, algorithm = "Lloyd")
      as.integer(km$cluster)
    },
    error = function(e) {
      km_fb <- kmeans(scale(X), centers = k, nstart = 25)
      as.integer(km_fb$cluster)
    }
  )
  
  list(cluster = res)
}



# -- 4. SIMULATION SETUP & EXECUTION --
# -------------------------------------
# Define parameters for the simulation.
radii_to_test <- seq(10, 1, by = -1) # Test from 10 down to 1
N_SHELLS <- 4
K_PER_SHELL <- 100
NOISE_SD <- 0.1
D_THRESHOLD <- 1.0 # Fixed distance threshold for the similarity graph

# Create an empty data frame to store results.
simulation_results <- tibble(
  max_radius = numeric(),
  estimated_k = integer()
)


# --- RUN THE SIMULATION ---
message("\nStarting simulation with Spectral Clustering. This may take a few minutes...")

for (R in radii_to_test) {
  message(paste0("Running: max_radius = ", R))
  
  # 1. Generate data (we only need the coordinates)
  shell_data <- generate_shell_clusters(N_SHELLS, K_PER_SHELL, R, NOISE_SD) %>%
    select(x, y, z)
  
  # 2. Run Gap Statistic with our custom spectral clustering function
  gap_stat <- clusGap(
    shell_data,
    FUNcluster = spectral_wrapper,
    K.max = N_SHELLS + 2, # Search up to k=6
    B = 50,
    d_threshold = D_THRESHOLD, # Pass the threshold to our wrapper
    verbose = FALSE
  )
  
  # 3. Estimate k and store the results
  estimated_k <- maxSE(f = gap_stat$Tab[, "gap"],
                       SE.f = gap_stat$Tab[, "SE.sim"],
                       method = "firstSEmax")
  
  simulation_results <- simulation_results %>%
    add_row(max_radius = R, estimated_k = estimated_k)
}

message("Simulation complete.")


# -- 6. VISUALIZE AND INTERPRET RESULTS --
# ----------------------------------------
# Plot the estimated number of clusters vs. the maximum radius.
final_plot <- ggplot(simulation_results, aes(x = max_radius, y = estimated_k)) +
  geom_hline(yintercept = N_SHELLS, color = "red", linetype = "dashed", linewidth = 1) +
  geom_line(color = "dodgerblue", linewidth = 1) +
  geom_point(color = "dodgerblue", size = 3) +
  scale_x_reverse(breaks = radii_to_test) +
  scale_y_continuous(breaks = 1:(N_SHELLS + 1)) +
  labs(
    title = "Spectral Clustering Performance on Concentric Shells",
    subtitle = "The red dashed line indicates the true number of clusters (4).",
    x = "Maximum Radius — Decreasing Shell Separation →",
    y = "Estimated Number of Clusters (k)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5)
  )

print(final_plot)

# ggsave("C:/Eneda_Research/FOXA1_ILC_RNA_DMSOE2/Spectral_performance.pdf",
#        plot = final_plot,
#        width = 12, height = 8, dpi = 300)
# 
# message("Plot saved")


# -- 7. INTERPRETATION OF RESULTS --
# ----------------------------------
# The plot shows the estimated number of clusters as the max_radius
# decreases. A smaller max_radius compresses the shells, making the
# distance between them smaller.
#
# Algorithm's Failure Point:
# The algorithm performs perfectly for max_radius values from 10 down
# to 8. In this range, the shells are sufficiently far apart. The average
# distance between adjacent shells is greater than our fixed `d_threshold` of 1.0.
# This allows the similarity graph to correctly connect points within a shell
# without incorrectly connecting points between shells.
#
# The failure point occurs around a max_radius of 7. At this point, the
# shells are so close that the distance between points on adjacent shells
# starts to fall below d_threshold. The similarity graph begins to form
# bridges between the shells, connecting them into a single large component.
# The Laplacian can no longer find a clean "cut" between them, and the Gap
# Statistic correctly identifies that the data now appears as one large cluster.
#
# - If d_threshold were smaller (e.g., 0.8): The algorithm would be
#   more "strict" about which points it considers "close". It would require
#   points to be nearer to form a connection. This would likely improve
#   performance, allowing the algorithm to distinguish shells at an even
#   smaller max_radius. 
#
# - If d_threshold were larger (e.g., 1.2): The algorithm would be
#   more "lenient", connecting points that are farther apart. This would
#   worsen performance. The algorithm would fail sooner, at a larger
#   max_radius. This is because it would start incorrectly connecting
#   adjacent shells much earlier in the process as they approach each other.

