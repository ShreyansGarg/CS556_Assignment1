#include <iostream>
#include <fstream>
#include <vector>
#include <queue>
#include <random>
#include <cmath>
#include <algorithm>
#include <iomanip>
#include <map>

// Define event types
enum EventType { ARRIVAL, DEPARTURE };

// Structure to represent an event
struct Event {
    double time;
    EventType type;
    
    // Constructor
    Event(double t, EventType typ) : time(t), type(typ) {}
    
    // Comparison operator for priority queue
    bool operator>(const Event& other) const {
        return time > other.time;
    }
};

// Class for the call center simulation
class CallCenterSimulation {
private:
    // Parameters
    double arrivalRate;                // Lambda (λ) - calls per hour
    double serviceRate;                // Mu (μ) - calls per hour
    int systemCapacity;                // c - including customer being served
    double simulationTime;             // in hours
    
    // Random number generators
    std::mt19937 rng;
    std::exponential_distribution<double> arrivalDistribution;
    std::exponential_distribution<double> serviceDistribution;
    
    // Simulation state
    double currentTime;
    int customersInSystem;
    int customersInQueue;
    std::priority_queue<Event, std::vector<Event>, std::greater<Event>> eventQueue;
    std::queue<double> serviceQueue;   // Queue of arrival times for waiting customers
    bool serverBusy;
    
    // Statistics
    int totalArrivals;
    int totalServices;
    int rejectedCustomers;
    double totalWaitingTime;
    double totalSystemTime;
    double totalBusyTime;
    double lastStateChangeTime;
    
    // For time-average calculations
    std::map<int, double> timeInState;  // Maps system state to cumulative time
    
    // Generate next arrival time
    double generateArrivalTime() {
        return arrivalDistribution(rng);
    }
    
    // Generate service time
    double generateServiceTime() {
        return serviceDistribution(rng);
    }
    
    // Update time spent in current state
    void updateTimeInState() {
        double timeSpent = currentTime - lastStateChangeTime;
        timeInState[customersInSystem] += timeSpent;
        lastStateChangeTime = currentTime;
    }
    
    // Process arrival event
    void processArrival() {
        // Schedule next arrival
        eventQueue.push(Event(currentTime + generateArrivalTime(), ARRIVAL));
        
        totalArrivals++;
        
        if (customersInSystem < systemCapacity) {
            // Customer can enter the system
            customersInSystem++;
            
            if (!serverBusy) {
                // Server is idle, begin service immediately
                serverBusy = true;
                // Schedule departure
                eventQueue.push(Event(currentTime + generateServiceTime(), DEPARTURE));
            } else {
                // Server is busy, join queue
                customersInQueue++;
                serviceQueue.push(currentTime);  // Record arrival time
            }
        } else {
            // System is full, customer is rejected
            rejectedCustomers++;
        }
    }
    
    // Process departure event
    void processDeparture() {
        // Customer completes service and leaves
        customersInSystem--;
        totalServices++;
        
        // Calculate time spent for this customer
        double serviceTime = 1.0 / serviceRate;  // Mean service time
        totalSystemTime += serviceTime;  // Add service time to system time
        
        // If there are waiting customers, start serving the next one
        if (customersInQueue > 0) {
            // Get arrival time of next customer
            double arrivalTime = serviceQueue.front();
            serviceQueue.pop();
            
            // Calculate waiting time (current time minus arrival time)
            double waitingTime = currentTime - arrivalTime;
            totalWaitingTime += waitingTime;
            totalSystemTime += waitingTime;  // Add waiting time to system time
            
            customersInQueue--;
            
            // Schedule next departure
            eventQueue.push(Event(currentTime + generateServiceTime(), DEPARTURE));
        } else {
            // No more customers in queue
            serverBusy = false;
        }
    }
    
    // Theoretical calculations for M/M/1/c
    void theoreticalCalculations(std::ofstream& outFile) {
        double rho = arrivalRate / serviceRate;
        
        // Calculate p0 (probability of empty system)
        double sum = 0;
        for (int i = 0; i <= systemCapacity; i++) {
            sum += std::pow(rho, i);
        }
        double p0 = 1.0 / sum;
        
        // Calculate probabilities for each state
        std::vector<double> p(systemCapacity + 1);
        for (int i = 0; i <= systemCapacity; i++) {
            p[i] = p0 * std::pow(rho, i);
        }
        
        // Performance metrics
        double rejectionProb = p[systemCapacity];
        double effectiveArrivalRate = arrivalRate * (1 - rejectionProb);
        double utilization = effectiveArrivalRate / serviceRate;
        
        // Average number in system
        double avgInSystem = 0;
        for (int i = 1; i <= systemCapacity; i++) {
            avgInSystem += i * p[i];
        }
        
        // Average number in queue
        double avgInQueue = 0;
        for (int i = 2; i <= systemCapacity; i++) {
            avgInQueue += (i - 1) * p[i];
        }
        
        // Average time in system and queue
        double avgTimeInSystem = avgInSystem / effectiveArrivalRate;
        double avgTimeInQueue = avgInQueue / effectiveArrivalRate;
        
        outFile << "\n---- Theoretical Results for M/M/1/" << systemCapacity << " Queue ----\n";
        outFile << "Theoretical Utilization: " << std::fixed << std::setprecision(4) << utilization * 100 << "%\n";
        outFile << "Theoretical Avg. Number in Queue: " << std::fixed << std::setprecision(4) << avgInQueue << std::endl;
        outFile << "Theoretical Avg. Number in System: " << std::fixed << std::setprecision(4) << avgInSystem << std::endl;
        outFile << "Theoretical Avg. Waiting Time (hours): " << std::fixed << std::setprecision(4) << avgTimeInQueue << " (" 
                << avgTimeInQueue * 60 << " minutes)\n";
        outFile << "Theoretical Avg. System Time (hours): " << std::fixed << std::setprecision(4) << avgTimeInSystem << " (" 
                << avgTimeInSystem * 60 << " minutes)\n";
        outFile << "Theoretical Rejection Probability: " << std::fixed << std::setprecision(4) << rejectionProb * 100 << "%\n";
        outFile << "Theoretical Probability of Full System: " << std::fixed << std::setprecision(4) << p[systemCapacity] * 100 << "%\n";
    }
    
public:
    // Constructor
    CallCenterSimulation(double lambda, double mu, int capacity, double simTime, unsigned seed = std::random_device{}())
        : arrivalRate(lambda),
          serviceRate(mu),
          systemCapacity(capacity),
          simulationTime(simTime),
          rng(seed),
          arrivalDistribution(lambda),
          serviceDistribution(mu),
          currentTime(0.0),
          customersInSystem(0),
          customersInQueue(0),
          serverBusy(false),
          totalArrivals(0),
          totalServices(0),
          rejectedCustomers(0),
          totalWaitingTime(0.0),
          totalSystemTime(0.0),
          totalBusyTime(0.0),
          lastStateChangeTime(0.0) {
        
        // Initialize state times to zero
        for (int i = 0; i <= systemCapacity; i++) {
            timeInState[i] = 0.0;
        }
        
        // Schedule first arrival
        eventQueue.push(Event(generateArrivalTime(), ARRIVAL));
    }
    
    // Run the simulation
    void runSimulation() {
        // Event loop
        while (!eventQueue.empty() && currentTime < simulationTime) {
            // Get next event
            Event event = eventQueue.top();
            eventQueue.pop();
            
            // Update time in state before processing event
            updateTimeInState();
            
            // Update current time
            currentTime = event.time;
            
            // Handle event based on type
            if (event.type == ARRIVAL) {
                processArrival();
            } else {  // DEPARTURE
                processDeparture();
            }
            
            // Update busy time
            if (serverBusy) {
                totalBusyTime += currentTime - lastStateChangeTime;
            }
        }
        
        // Add final state time
        updateTimeInState();
    }
    
    // Calculate performance metrics
    struct PerformanceMetrics {
        double avgWaitingTime;       // Average waiting time in hours
        double avgSystemTime;        // Average time in system in hours
        double utilization;          // Agent utilization rate
        double avgQueueLength;       // Average number of customers in queue
        double fullSystemProb;       // Probability of system being full
        double rejectionProb;        // Probability of customer rejection
    };
    
    PerformanceMetrics calculateMetrics() {
        PerformanceMetrics metrics;
        
        // Average waiting time in hours
        metrics.avgWaitingTime = totalServices > 0 ? totalWaitingTime / totalServices : 0;
        
        // Average time in system in hours (waiting + service)
        metrics.avgSystemTime = totalServices > 0 ? totalSystemTime / totalServices : 0;
        
        // Agent utilization rate
        metrics.utilization = totalBusyTime / currentTime;
        
        // Average queue length (time-weighted)
        metrics.avgQueueLength = 0;
        for (int i = 2; i <= systemCapacity; i++) {
            metrics.avgQueueLength += (i - 1) * (timeInState[i] / currentTime);
        }
        
        // Probability of system being full
        metrics.fullSystemProb = timeInState[systemCapacity] / currentTime;
        
        // Probability of customer rejection
        metrics.rejectionProb = (totalArrivals + rejectedCustomers) > 0 ? 
                                static_cast<double>(rejectedCustomers) / (totalArrivals + rejectedCustomers) : 0;
        
        return metrics;
    }
    
    // Analyze system with varying capacities
    void analyzeVaryingCapacity(std::ofstream& outFile) {
        outFile << "\n---- Analysis of Varying System Capacity ----\n";
        outFile << std::setw(10) << "Capacity" 
                << std::setw(15) << "Avg Wait(min)" 
                << std::setw(15) << "Utilization(%)" 
                << std::setw(15) << "Rejection(%)" << std::endl;
        
        for (int capacity = 3; capacity <= 7; capacity++) {
            CallCenterSimulation sim(arrivalRate, serviceRate, capacity, simulationTime);
            sim.runSimulation();
            auto metrics = sim.calculateMetrics();
            
            outFile << std::setw(10) << capacity
                    << std::setw(15) << std::fixed << std::setprecision(2) << metrics.avgWaitingTime * 60 // convert to minutes
                    << std::setw(15) << std::fixed << std::setprecision(2) << metrics.utilization * 100 // as percentage
                    << std::setw(15) << std::fixed << std::setprecision(2) << metrics.rejectionProb * 100 // as percentage
                    << std::endl;
            
            // Create data file for plotting
            std::ofstream plotData("capacity_" + std::to_string(capacity) + "_data.txt");
            plotData << "Avg_Waiting_Time_Minutes " << metrics.avgWaitingTime * 60 << std::endl;
            plotData << "Utilization " << metrics.utilization * 100 << std::endl;
            plotData << "Rejection_Probability " << metrics.rejectionProb * 100 << std::endl;
            plotData.close();
        }
    }
    
    // Generate report with all results
    void generateReport(std::ofstream& outFile) {
        // Calculate metrics
        auto metrics = calculateMetrics();
        
        // Write simulation parameters
        outFile << "==== Call Center M/M/1/" << systemCapacity << " Queue Simulation Report ====\n\n";
        outFile << "Simulation Parameters:\n";
        outFile << "- Arrival Rate (λ): " << arrivalRate << " calls/hour\n";
        outFile << "- Service Rate (μ): " << serviceRate << " calls/hour\n";
        outFile << "- System Capacity (c): " << systemCapacity << " (including customer being served)\n";
        outFile << "- Simulation Time: " << simulationTime << " hours\n\n";
        
        // Write performance metrics
        outFile << "Performance Metrics:\n";
        outFile << "- Average Waiting Time in Queue: " << std::fixed << std::setprecision(4) << metrics.avgWaitingTime << " hours (" 
                << metrics.avgWaitingTime * 60 << " minutes)\n";
        outFile << "- Average Time in System: " << std::fixed << std::setprecision(4) << metrics.avgSystemTime << " hours (" 
                << metrics.avgSystemTime * 60 << " minutes)\n";
        outFile << "- Agent Utilization Rate: " << std::fixed << std::setprecision(2) << metrics.utilization * 100 << "%\n";
        outFile << "- Average Number of Customers in Queue: " << std::fixed << std::setprecision(4) << metrics.avgQueueLength << "\n";
        outFile << "- Probability of System Being Full: " << std::fixed << std::setprecision(2) << metrics.fullSystemProb * 100 << "%\n";
        outFile << "- Probability of Customer Rejection: " << std::fixed << std::setprecision(2) << metrics.rejectionProb * 100 << "%\n\n";
        
        outFile << "Simulation Statistics:\n";
        outFile << "- Total Customers Arrived: " << totalArrivals + rejectedCustomers << "\n";
        outFile << "- Customers Served: " << totalServices << "\n";
        outFile << "- Customers Rejected: " << rejectedCustomers << "\n\n";
        
        // Calculate and write theoretical results
        theoreticalCalculations(outFile);
        
        // Write analysis of varying capacity
        analyzeVaryingCapacity(outFile);
        
        // Write recommendations
        outFile << "\n==== Recommendations ====\n\n";
        
        // Determine optimal capacity for < 5 min waiting time
        bool foundOptimal = false;
        int optimalCapacity = 0;
        double minRejectionRate = 100.0;
        
        for (int capacity = 3; capacity <= 7; capacity++) {
            CallCenterSimulation sim(arrivalRate, serviceRate, capacity, simulationTime);
            sim.runSimulation();
            auto capacityMetrics = sim.calculateMetrics();
            
            // Check if waiting time is less than 5 minutes
            if (capacityMetrics.avgWaitingTime * 60 < 5.0) {
                if (!foundOptimal || capacityMetrics.rejectionProb < minRejectionRate) {
                    optimalCapacity = capacity;
                    minRejectionRate = capacityMetrics.rejectionProb;
                    foundOptimal = true;
                }
            }
        }
        
        outFile << "1. Optimal System Capacity:\n";
        if (foundOptimal) {
            outFile << "   Based on the simulation results, a system capacity of " << optimalCapacity 
                    << " is recommended to maintain an average waiting time of less than 5 minutes "
                    << "while minimizing rejections.\n\n";
        } else {
            outFile << "   None of the tested capacities (3-7) achieved an average waiting time of less than 5 minutes. "
                    << "Consider increasing service rate or system capacity beyond 7.\n\n";
        }
        
        outFile << "2. Capacity Planning Strategies:\n";
        outFile << "   - Consider implementing a dynamic capacity system that can adjust based on demand patterns.\n";
        outFile << "   - Monitor the trade-off between waiting time and rejection probability regularly.\n";
        outFile << "   - For periods of high call volume, consider temporary increases in system capacity.\n\n";
        
        outFile << "3. Performance Optimization:\n";
        outFile << "   - If the agent utilization is consistently high (>85%), consider adding another agent.\n";
        outFile << "   - Implement a callback system for customers when the system is near capacity.\n";
        outFile << "   - Consider specialized training for the agent to reduce service time without compromising quality.\n";
        outFile << "   - Implement a pre-recorded information system to handle common queries before connecting to the agent.\n\n";
        
        outFile << "4. Assumptions:\n";
        outFile << "   - Poisson arrival process (exponential inter-arrival times)\n";
        outFile << "   - Exponential service times\n";
        outFile << "   - First-come, first-served queue discipline\n";
        outFile << "   - Stationary arrival and service rates (no time-of-day variations)\n";
        outFile << "   - No customer abandonment (reneging or balking)\n";
        outFile << "   - Single service channel\n\n";
        
        outFile << "==== End of Report ====\n";
    }
};

int main() {
    // Parameters
    double arrivalRate = 20.0;      // λ (lambda) - calls per hour
    double serviceRate = 24.0;      // μ (mu) - calls per hour
    int systemCapacity = 5;         // c - including customer being served
    double simulationTime = 1000.0; // hours (long enough for steady-state)
    
    // Create output file
    std::ofstream outFile("call_center_simulation_results.txt");
    if (!outFile.is_open()) {
        std::cerr << "Error: Unable to open output file." << std::endl;
        return 1;
    }
    
    // Create and run simulation
    CallCenterSimulation simulation(arrivalRate, serviceRate, systemCapacity, simulationTime);
    simulation.runSimulation();
    
    // Generate report
    simulation.generateReport(outFile);
    
    // Close file
    outFile.close();
    
    // Print completion message to terminal
    std::cout << "Simulation completed. Results have been written to 'call_center_simulation_results.txt'" << std::endl;
    
    return 0;
}