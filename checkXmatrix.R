# Generate output for X1_data
y1 <- calculate_NLES5_for_sobol(X1_data)

# Generate output for X2_data (soboljansen also uses this)
y2 <- calculate_NLES5_for_sobol(X2_data)

# --- Now, let's check y1 and y2 ---

# Check for any NaN, Inf, or NA values:
print(paste("Any NaN in y1:", any(is.nan(y1))))
print(paste("Any Inf in y1:", any(is.infinite(y1))))
print(paste("Any NA in y1:", any(is.na(y1))))
print(summary(y1))

print(paste("Any NaN in y2:", any(is.nan(y2))))
print(paste("Any Inf in y2:", any(is.infinite(y2))))
print(paste("Any NA in y2:", any(is.na(y2))))
print(summary(y2))

# Check the variance of the output. If it's zero, indices will be NaN.
print(paste("Variance of y1:", var(y1, na.rm = TRUE)))
print(paste("Variance of y2:", var(y2, na.rm = TRUE)))

# Check how many unique values are in the output. If it's 1, the output is constant.
print(paste("Number of unique values in y1:", length(unique(y1))))
print(paste("Number of unique values in y2:", length(unique(y2))))

# View the first few and last few values
print(head(y1))
print(tail(y1))


# Check summaries for NAs or strange values
summary(X1_data)

print("Summary of X2_data:")
summary(X2_data)

# Check if any column in X1_data or X2_data has zero variance
print(sapply(X1_data, var, na.rm = TRUE))

# Check for NAs in the sampled data (though na.omit should handle this in sampling)
print(paste("Any NA in X1_data:", any(is.na(X1_data))))
print(paste("Any NA in X2_data:", any(is.na(X2_data))))

