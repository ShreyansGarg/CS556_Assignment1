import matplotlib.pyplot as plt
import numpy as np
import os

# Function to read data from files
def read_metrics_data(capacity_range):
    waiting_times = []
    utilizations = []
    rejection_rates = []
    
    for capacity in capacity_range:
        filename = f"capacity_{capacity}_data.txt"
        if os.path.exists(filename):
            with open(filename, 'r') as file:
                data = {}
                for line in file:
                    key, value = line.strip().split(' ')
                    data[key] = float(value)
                    
            waiting_times.append(data.get('Avg_Waiting_Time_Minutes', 0))
            utilizations.append(data.get('Utilization', 0))
            rejection_rates.append(data.get('Rejection_Probability', 0))
        else:
            print(f"Warning: File {filename} not found. Using 0 values.")
            waiting_times.append(0)
            utilizations.append(0)
            rejection_rates.append(0)
            
    return waiting_times, utilizations, rejection_rates

def plot_metrics(capacity_range):
    waiting_times, utilizations, rejection_rates = read_metrics_data(capacity_range)
    
    # Create figure with 3 subplots
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(10, 12))
    fig.suptitle('Call Center Performance Metrics vs. System Capacity', fontsize=16)
    
    # Plot 1: Average Waiting Time
    ax1.plot(capacity_range, waiting_times, 'o-', color='blue', linewidth=2, markersize=8)
    ax1.set_ylabel('Average Waiting Time (minutes)')
    ax1.set_xlabel('System Capacity')
    ax1.grid(True)
    ax1.set_title('Average Waiting Time vs. System Capacity')
    
    # Plot 2: Agent Utilization
    ax2.plot(capacity_range, utilizations, 'o-', color='green', linewidth=2, markersize=8)
    ax2.set_ylabel('Agent Utilization (%)')
    ax2.set_xlabel('System Capacity')
    ax2.grid(True)
    ax2.set_title('Agent Utilization vs. System Capacity')
    
    # Plot 3: Rejection Probability
    ax3.plot(capacity_range, rejection_rates, 'o-', color='red', linewidth=2, markersize=8)
    ax3.set_ylabel('Rejection Probability (%)')
    ax3.set_xlabel('System Capacity')
    ax3.grid(True)
    ax3.set_title('Rejection Probability vs. System Capacity')
    
    # Adjust layout
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    
    # Save plot
    plt.savefig('call_center_metrics.png', dpi=300)
    
    # Create combined plot
    plt.figure(figsize=(12, 8))
    plt.plot(capacity_range, waiting_times, 'o-', color='blue', linewidth=2, label='Avg. Waiting Time (min)')
    
    # Create a second y-axis for percentage metrics
    ax4 = plt.gca()
    ax5 = ax4.twinx()
    
    ax5.plot(capacity_range, utilizations, 's-', color='green', linewidth=2, label='Utilization (%)')
    ax5.plot(capacity_range, rejection_rates, '^-', color='red', linewidth=2, label='Rejection Rate (%)')
    
    ax4.set_xlabel('System Capacity', fontsize=12)
    ax4.set_ylabel('Average Waiting Time (minutes)', color='blue', fontsize=12)
    ax5.set_ylabel('Percentage (%)', color='black', fontsize=12)
    
    # Add a legend
    lines1, labels1 = ax4.get_legend_handles_labels()
    lines2, labels2 = ax5.get_legend_handles_labels()
    ax5.legend(lines1 + lines2, labels1 + labels2, loc='upper right')
    
    plt.title('Call Center Performance Metrics vs. System Capacity', fontsize=14)
    plt.grid(True)
    plt.tight_layout()
    
    # Save combined plot
    plt.savefig('call_center_combined_metrics.png', dpi=300)
    
    print("Plots generated: call_center_metrics.png and call_center_combined_metrics.png")

if __name__ == "__main__":
    # Define range of capacities to plot
    capacity_range = range(3, 8)  # Capacities 3 to 7
    
    plot_metrics(capacity_range)
