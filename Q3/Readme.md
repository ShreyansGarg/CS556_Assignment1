# Coffee Shop Queue Simulation (M/M/1 & M/M/c with Reneging)

## Overview

This project simulates a **single-server (M/M/1)** and **multi-server (M/M/c)** queue model for a coffee shop, where customers arrive randomly and are served on a **first-come, first-served (FCFS)** basis. The simulation also includes an analysis of customer reneging (abandonment) when wait times exceed a threshold.

## Problem Statement

A small coffee shop operates with a **single barista** serving customers. The objective is to analyze the system’s performance using an **M/M/1 queue model** and evaluate:

- **Average waiting time**
- **Service efficiency**
- **Customer abandonment (reneging)**
- **Impact of adding a second barista (M/M/2 model)**
- **Revenue lost due to abandoned customers**

## System Description

1. **Customers arrive randomly**, following a **Poisson process** with arrival rate **λ = 10 customers/hour**.
2. **Service times follow an exponential distribution** with service rate **μ = 15 customers/hour**.
3. Customers who arrive when the barista is busy **wait in a queue**.
4. If a customer’s wait time exceeds **5 minutes (reneging threshold)**, they **abandon** the queue.
5. The simulation runs until **500 customers** have been processed.

## Files in This Project

### 1. `q3reneging.cpp` (Simulation Code)

- Implements a **discrete event simulation (DES)** for M/M/1 and M/M/c queue systems.
- Models **customer arrivals, service, and departures**.
- Simulates **reneging** when wait times exceed **5 minutes**.
- Compares performance under different conditions:
  - **Single barista (M/M/1)**
  - **Faster service (μ = 20 customers/hour)**
  - **Slower service (μ = 12 customers/hour)**
  - **Two baristas (M/M/2)**
  - **Customer abandonment impact**
- Outputs results to `coffee_shop_simulation_results.txt`.

### 2. `Makefile` (Build & Run Automation)

This Makefile automates:

- **Compilation** of `q3reneging.cpp`
- **Execution** of the simulation

#### Usage:

```bash
make        # Compile and run the simulation
make clean  # Remove compiled files
```

---

## How to Run the Simulation

### **Step 1: Compile & Run**

```bash
make
```

This compiles `q3reneging.cpp` and runs the simulation, generating output files with performance metrics.

### **Step 2: View the Results**

- **Simulation Report (Without Reneging):** `coffee_shop_simulation_results.txt`
- **Simulation Report (With Reneging):** `coffee_shop_simulation_with_fixed_reneging_results.txt`

---

## Performance Metrics Computed

1. **Average Waiting Time in Queue** – Time spent waiting before ordering.
2. **Average Time in System** – Total time in shop (waiting + service).
3. **Utilization Factor** – Percentage of time the barista is busy.
4. **Idle Time of the Barista** – Time the barista is not serving.
5. **Maximum Queue Length** – Longest waiting line observed.
6. **Probability of an Empty Queue** – Percentage of time with no waiting customers.
7. **Peak Hour Analysis** – Identifies the busiest period.
8. **Customer Reneging Rate** – Percentage of customers who leave due to long wait times.
9. **Revenue Lost Due to Reneging** – Missed earnings from abandoned customers.

---

## Business Insights & Recommendations

1. **Speeding up service (μ = 20 customers/hour) significantly reduces wait times and abandonment.**
2. **Adding a second barista (M/M/2) improves efficiency and prevents long queues.**
3. **Customer abandonment leads to revenue losses** – strategies like mobile ordering can help.
4. **Dynamic staffing (hiring extra baristas during peak hours) can reduce wait times.**
5. **Implementing a loyalty program could incentivize customers to stay even with longer wait times.**

---

## Assumptions

- Customers **arrive randomly (Poisson process)**.
- Service times **follow an exponential distribution**.
- Customers are **served in a first-come, first-served (FCFS) manner**.
- If **reneging is enabled, customers abandon after waiting 5 minutes**.
- No customer retries – abandoned customers do not return.

---

## Conclusion

This simulation helps understand the **impact of service speed, customer reneging, and additional staff** in a coffee shop. By analyzing key metrics, the shop can make **data-driven decisions to improve service quality, reduce lost revenue, and optimize operations**.

