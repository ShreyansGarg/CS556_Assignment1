# Bank Teller M/M/c Queue Simulation

## Overview
This project simulates a bank with a varying number of tellers throughout the workday using an **M/M/c queue model**. Customers arrive based on a Poisson process, and service times follow an exponential distribution. The goal is to evaluate the system's performance and optimize teller allocation.

## Simulation Details
- **Customer arrival rate:** 40 customers per hour (Poisson process)
- **Service time:** Exponentially distributed, averaging 4 minutes per customer
- **Time Periods & Tellers:**
  - 9:00 AM - 11:00 AM: **2 tellers**
  - 11:00 AM - 2:00 PM: **4 tellers** (peak hours)
  - 2:00 PM - 5:00 PM: **3 tellers**
- **Total simulation duration:** 8 hours (9:00 AM - 5:00 PM)

## Performance Metrics
The simulation calculates:
1. **Average waiting time in queue**
2. **Average time spent in the system (waiting + service)**
3. **Teller utilization rate**
4. **Average number of customers in queue**
5. **Probability of all tellers being busy**
6. **Queue length and waiting time plots**

## File Structure
- **`q1.cpp`**: C++ program simulating the queue system
- **`queue-waiting-time-plots.py`**: Python script for visualizing queue data
- **`Makefile`**: Automates the build and run process
- **`report.txt`**: Outputs the simulation results
- **Output CSVs**: `queue_length.csv`, `waiting_time.csv`
- **Plots**:
  - `queue_and_waiting_time_plots.png`
  - `queue_length_vs_waiting_time.png`
  - `avg_waiting_time_by_interval.png`

## Running the Simulation
### **Using Makefile**
To compile and run the program, execute:
```sh
make
```
To clean compiled files:
```sh
make clean
```

### **Manual Execution**
1. Compile and run the C++ program:
   ```sh
   g++ q1.cpp -o q1 && ./q1
   ```
   Results are stored in `report.txt`.
2. Generate visualizations using Python:
   ```sh
   python3 queue-waiting-time-plots.py
   ```
   This will produce queue length and waiting time plots.

## Findings and Recommendations
- **Peak hours (11 AM - 2 PM)** show the highest queue length and waiting times.
- Increasing tellers in high-wait periods can improve service efficiency.
- Optimal teller allocation can reduce customer waiting times while maintaining efficient utilization.




---


