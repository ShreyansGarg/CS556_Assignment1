# CS556_Assignment1

# Combined Queue Simulation Project (M/M/1 & M/M/c)

## Overview

This project simulates queue management for different service-based environments using **M/M/1** and **M/M/c** and **M/M/1/c** queue models. The simulations analyze customer arrivals, service efficiency, waiting times, and system performance under varying conditions. Three primary scenarios are covered:

1. **Bank Teller Queue (M/M/c Model)** – Varying teller allocation throughout the workday.
2. **Call Center Queue (M/M/1/c Model)** – Limited capacity with customer rejections.
3. **Coffee Shop Queue (M/M/1 & M/M/c with Reneging)** – Customer abandonment due to long wait times.

## System Descriptions

### 1. Bank Teller Queue (M/M/c)
- **Customer arrival rate:** 40 customers per hour (Poisson process).
- **Service time:** Exponentially distributed, averaging 4 minutes per customer.
- **Teller allocation:**
  - **9:00 AM - 11:00 AM:** 2 tellers.
  - **11:00 AM - 2:00 PM:** 4 tellers (peak hours).
  - **2:00 PM - 5:00 PM:** 3 tellers.
- **Simulation duration:** 8 hours (9:00 AM - 5:00 PM).

### 2. Call Center Queue (M/M/1/c)
- **Arrival rate (λ):** 20 calls per hour.
- **Service rate (μ):** 24 calls per hour (single agent).
- **System capacity (c):** 5 (including the customer being served).
- **Performance metrics include:**
  - Average waiting time.
  - Probability of system being full.
  - Customer rejection rate.
  - Utilization of the agent.
- **Capacity variation:** The system is tested with capacities from **3 to 7**.

### 3. Coffee Shop Queue (M/M/1 & M/M/c with Reneging)
- **Customer arrival rate:** 10 customers per hour.
- **Service rate:** 15 customers per hour.
- **Queue discipline:** First-come, first-served (FCFS).
- **Customer abandonment:** If waiting time exceeds 5 minutes, customers leave.
- **Comparison of different scenarios:**
  - Single barista (M/M/1).
  - Faster service (μ = 20 customers/hour).
  - Slower service (μ = 12 customers/hour).
  - Two baristas (M/M/2).
  - Revenue impact due to abandoned customers.

---

## Files in This Project

### 1. **C++ Simulation Programs**
- `q1.cpp` (Bank Teller Queue Simulation).
- `q2.cpp` (Call Center Queue Simulation).
- `q3reneging.cpp` (Coffee Shop Queue with Reneging).

### 2. **Python Visualization Scripts**
- `queue-waiting-time-plots.py` (Bank Teller queue visualization).
- `plot-call-center-metrics.py` (Call Center queue visualization).

### 3. **Makefile (Build & Run Automation)**
Automates:
- **Compilation** of all C++ programs.
- **Execution** of simulations.
- **Data visualization using Python.**
- 
### The methods for using makefile for running the code are in the individual Readme of each problem.


### ** View the Results**
- **Bank Teller Results:** `report.txt`.
- **Call Center Results:** `call_center_simulation_results.txt`.
- **Coffee Shop Results:** `coffee_shop_simulation_results.txt`.
- **Visualization Outputs:**
  - `queue_and_waiting_time_plots.png`
  - `call_center_metrics.png`
  - `coffee_shop_metrics.png`

---

## Performance Metrics Computed

1. **Average Waiting Time in Queue** – How long customers wait before being served.
2. **Average Time in System** – Total time in the system (waiting + service).
3. **Utilization Factor** – Percentage of time servers are busy.
4. **Idle Time** – Percentage of time servers are not serving customers.
5. **Maximum Queue Length** – The longest observed queue length.
6. **Probability of System Being Full (Call Center)** – Chance that all slots are occupied.
7. **Customer Reneging Rate (Coffee Shop)** – Percentage of customers who leave.
8. **Revenue Lost Due to Abandonment (Coffee Shop)** – Missed earnings from customers who left due to long wait times.

---

## Business Insights & Recommendations

1. **Dynamic Teller Scheduling** – Increase tellers during peak hours to reduce waiting times.
2. **Call Center Capacity Planning** – Adjust system capacity to balance customer rejection and agent utilization.
3. **Coffee Shop Optimization** – Add a second barista or implement mobile ordering to minimize lost revenue.
4. **Data-Driven Decision Making** – Use real-time queue monitoring to adjust resources dynamically.

---

## Assumptions

- Customers **arrive randomly (Poisson process)**.
- Service times **follow an exponential distribution**.
- **First-come, first-served (FCFS)** queuing discipline.
- No priority customers.
- No service interruptions or external influences.

---

## Conclusion

This project provides insights into **service optimization in banks, call centers, and coffee shops** using **M/M/1** and **M/M/c** queue models. The simulations help businesses make **data-driven decisions to enhance efficiency, reduce wait times, and improve customer experience.**

