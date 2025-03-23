import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Load the data
queue_length_df = pd.read_csv('queue_length.csv')
waiting_time_df = pd.read_csv('waiting_time.csv')

# Create a figure with two subplots
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10))

# Plot 1: Queue Length Over Time
ax1.plot(queue_length_df['Time'], queue_length_df['Value'], 'b-', linewidth=1)
ax1.set_title('Queue Length Over Time', fontsize=14)
ax1.set_xlabel('Time (hours)', fontsize=12)
ax1.set_ylabel('Queue Length', fontsize=12)
ax1.grid(True, linestyle='--', alpha=0.7)

# Add some statistics to the first plot
max_queue = queue_length_df['Value'].max()
avg_queue = queue_length_df['Value'].mean()
ax1.axhline(y=avg_queue, color='r', linestyle='--', alpha=0.8)
ax1.text(0.02, max_queue * 0.95, f'Maximum Queue Length: {max_queue}', 
         fontsize=10, bbox=dict(facecolor='white', alpha=0.8))
ax1.text(0.02, max_queue * 0.85, f'Average Queue Length: {avg_queue:.2f}', 
         fontsize=10, bbox=dict(facecolor='white', alpha=0.8))

# Plot 2: Waiting Time Over Time
ax2.plot(waiting_time_df['Time'], waiting_time_df['Value'], 'g-', linewidth=1)
ax2.set_title('Waiting Time Over Time', fontsize=14)
ax2.set_xlabel('Time (hours)', fontsize=12)
ax2.set_ylabel('Waiting Time (hours)', fontsize=12)
ax2.grid(True, linestyle='--', alpha=0.7)

# Add some statistics to the second plot
max_wait = waiting_time_df['Value'].max()
avg_wait = waiting_time_df['Value'].mean()
non_zero_wait = waiting_time_df[waiting_time_df['Value'] > 0]['Value']
avg_non_zero_wait = non_zero_wait.mean() if len(non_zero_wait) > 0 else 0

ax2.axhline(y=avg_wait, color='r', linestyle='--', alpha=0.8)
ax2.text(0.02, max_wait * 0.95, f'Maximum Waiting Time: {max_wait:.4f} hours', 
         fontsize=10, bbox=dict(facecolor='white', alpha=0.8))
ax2.text(0.02, max_wait * 0.85, f'Average Waiting Time: {avg_wait:.4f} hours', 
         fontsize=10, bbox=dict(facecolor='white', alpha=0.8))
ax2.text(0.02, max_wait * 0.75, f'Average Non-Zero Waiting Time: {avg_non_zero_wait:.4f} hours', 
         fontsize=10, bbox=dict(facecolor='white', alpha=0.8))

# Adjust layout
plt.tight_layout()

# Save the plots
plt.savefig('queue_and_waiting_time_plots.png', dpi=300)

# Show the plots
plt.show()

# Create a scatter plot of queue length vs waiting time
# First, create a combined dataframe with all unique timestamps
all_times = pd.concat([queue_length_df['Time'], waiting_time_df['Time']]).unique()
all_times.sort()
combined_df = pd.DataFrame({'Time': all_times})

# Interpolate queue length and waiting time for all timestamps
combined_df['Queue_Length'] = np.interp(combined_df['Time'], 
                                       queue_length_df['Time'], 
                                       queue_length_df['Value'])

combined_df['Waiting_Time'] = np.interp(combined_df['Time'], 
                                      waiting_time_df['Time'], 
                                      waiting_time_df['Value'])

# Create a scatter plot of queue length vs waiting time
plt.figure(figsize=(10, 6))
plt.scatter(combined_df['Queue_Length'], combined_df['Waiting_Time'], alpha=0.6)
plt.title('Queue Length vs Waiting Time', fontsize=14)
plt.xlabel('Queue Length', fontsize=12)
plt.ylabel('Waiting Time (hours)', fontsize=12)
plt.grid(True, linestyle='--', alpha=0.7)

# Add a trend line
mask = combined_df['Queue_Length'] > 0  # Only use non-zero queue length
if sum(mask) > 1:  # Check if we have enough data for a trend line
    z = np.polyfit(combined_df.loc[mask, 'Queue_Length'], 
                 combined_df.loc[mask, 'Waiting_Time'], 1)
    p = np.poly1d(z)
    plt.plot(np.sort(combined_df.loc[mask, 'Queue_Length']), 
           p(np.sort(combined_df.loc[mask, 'Queue_Length'])), 
           "r--", linewidth=2)
    plt.text(0.5, max(combined_df['Waiting_Time']) * 0.9, 
           f'Trend: y = {z[0]:.4f}x + {z[1]:.4f}', 
           fontsize=10, bbox=dict(facecolor='white', alpha=0.8))

plt.tight_layout()
plt.savefig('queue_length_vs_waiting_time.png', dpi=300)
plt.show()

# Create a bar graph for the average waiting time in 3-hour intervals
# Define time intervals
intervals = [(0, 2), (2, 5), (5, 8)]
interval_labels = ['0-2 hours', '2-5 hours', '5-8 hours']

# Calculate average waiting time for each interval
avg_waiting_times = []
for start, end in intervals:
    mask = (waiting_time_df['Time'] >= start) & (waiting_time_df['Time'] < end)
    interval_data = waiting_time_df.loc[mask, 'Value']
    avg_waiting_times.append(interval_data.mean())

# Create a bar graph
plt.figure(figsize=(10, 6))
bars = plt.bar(interval_labels, avg_waiting_times, color=['blue', 'green', 'orange'])

# Add values on top of the bars
for bar in bars:
    height = bar.get_height()
    plt.text(bar.get_x() + bar.get_width()/2., height + 0.001,
             f'{height:.4f}',
             ha='center', va='bottom', fontsize=10)

plt.title('Average Waiting Time by Time Interval', fontsize=14)
plt.xlabel('Time Interval', fontsize=12)
plt.ylabel('Average Waiting Time (hours)', fontsize=12)
plt.grid(True, linestyle='--', alpha=0.7, axis='y')
plt.tight_layout()
plt.savefig('avg_waiting_time_by_interval.png', dpi=300)
plt.show()