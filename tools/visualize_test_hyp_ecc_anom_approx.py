import pandas as pd
import matplotlib.pyplot as plt

# Step 1: Load the CSV file
data = pd.read_csv('./out/test_hyp_ecc_anom_approx.csv')

# Step 2: Group by 'orbit'
grouped_data = data.groupby('orbit')

# Step 3: Plot for each orbit
for type_name, group in grouped_data:
    plt.figure()  # Create a new figure for each orbit
    plt.title(f'Graph for {type_name}')
    plt.xlabel('iter')
    plt.ylabel('Values')

    # Plot the 'approx', 'real', and 'error' columns
    plt.plot(group['iter'], group['approx'], label='approx')
    plt.plot(group['iter'], group['real'], label='real')
    # plt.plot(group['iter'], group['error'], label='error')

    # Add a legend and show the plot
    plt.legend()
    plt.grid(True)
    plt.show()
