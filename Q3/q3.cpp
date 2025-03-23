#include <iostream>
#include <fstream>
#include <queue>
#include <random>
#include <cmath>
#include <iomanip>
#include <vector>
#include<map>
#include <algorithm>

struct Customer {
    double arrivalTime;
    double serviceStartTime;
    double departureTime;
    double maximumWaitTime;  // For reneging simulation
    double revenue;
    bool abandoned;
};

struct Event {
    enum Type { ARRIVAL, DEPARTURE, ABANDONMENT };
    
    Type type;
    double time;
    int customerId;
    int serverId;
    
    // For priority queue ordering (earliest events first)
    bool operator>(const Event& other) const {
        return time > other.time;
    }
};

struct ServerStats {
    double utilization;
    double idleTime;
    double avgWaitingTime;
    double avgSystemTime;
    int maxQueueLength;
    double emptyQueueProbability;
    double totalRevenue;
    double peakHour;
    int peakQueueLength;
    int abandonedCustomers;
    double lostRevenue;
};

// Function to generate random numbers from an exponential distribution
double generateExponential(double rate, std::mt19937& gen) {
    std::exponential_distribution<double> distribution(rate);
    return distribution(gen);
}

// Function to simulate queue with specified parameters and reneging behavior
ServerStats simulateQueue(double lambda, double mu, int numServers, int numCustomers, 
                         bool enableReneging, double patienceScale, 
                         bool printDetails, std::ofstream& outFile, std::mt19937& gen) {
    // Simulation variables
    std::vector<Customer> customers;
    std::vector<int> waitingQueue;
    std::vector<bool> serverBusy(numServers, false);
    std::priority_queue<Event, std::vector<Event>, std::greater<Event>> eventQueue;
    
    double currentTime = 0.0;
    double totalBusyTime = 0.0;
    int maxQueueLength = 0;
    int emptyQueueCount = 0;
    double totalRevenue = 0.0;
    double lostRevenue = 0.0;
    int completedCustomers = 0;
    int abandonedCustomers = 0;
    
    // For peak hour analysis
    std::map<int, int> hourlyQueueCounts;
    int peakQueueLength = 0;
    double peakHour = 0.0;
    
    // Schedule first arrival
    Event firstArrival;
    firstArrival.type = Event::ARRIVAL;
    firstArrival.time = generateExponential(lambda, gen);
    firstArrival.customerId = 0;
    eventQueue.push(firstArrival);
    
    // Process events until required number of customers have been served or abandoned
    while (completedCustomers + abandonedCustomers < numCustomers && !eventQueue.empty()) {
        // Get next event
        Event event = eventQueue.top();
        eventQueue.pop();
        
        // Update simulation time
        double timeDelta = event.time - currentTime;
        currentTime = event.time;
        
        // Update busy time for servers
        for (int i = 0; i < numServers; i++) {
            if (serverBusy[i]) {
                totalBusyTime += timeDelta;
            }
        }
        
        // Record queue stats
        int hourBucket = static_cast<int>(currentTime);
        hourlyQueueCounts[hourBucket] = std::max(hourlyQueueCounts[hourBucket], static_cast<int>(waitingQueue.size()));
        
        // Update empty queue count
        if (waitingQueue.empty()) {
            emptyQueueCount++;
        }
        
        // Process the event
        if (event.type == Event::ARRIVAL) {
            // Create new customer
            Customer newCustomer;
            newCustomer.arrivalTime = currentTime;
            newCustomer.revenue = 5.0;  // $5 per customer
            newCustomer.abandoned = false;
            
            if (enableReneging) {
                // Set maximum wait time for possible abandonment
                newCustomer.maximumWaitTime = generateExponential(1.0 / patienceScale, gen);
            } else {
                newCustomer.maximumWaitTime = std::numeric_limits<double>::infinity();
            }
            
            // Add customer to the system
            customers.push_back(newCustomer);
            int customerId = customers.size() - 1;
            
            // Check if a server is available
            bool foundServer = false;
            for (int i = 0; i < numServers; i++) {
                if (!serverBusy[i]) {
                    // Server available, start service
                    serverBusy[i] = true;
                    customers[customerId].serviceStartTime = currentTime;
                    double serviceTime = generateExponential(mu, gen);
                    customers[customerId].departureTime = currentTime + serviceTime;
                    
                    // Schedule departure event
                    Event departure;
                    departure.type = Event::DEPARTURE;
                    departure.time = customers[customerId].departureTime;
                    departure.customerId = customerId;
                    departure.serverId = i;
                    eventQueue.push(departure);
                    
                    foundServer = true;
                    break;
                }
            }
            
            if (!foundServer) {
                // All servers busy, join the waiting queue
                waitingQueue.push_back(customerId);
                
                // Update maximum queue length
                if (waitingQueue.size() > maxQueueLength) {
                    maxQueueLength = waitingQueue.size();
                    
                    // Update peak hour information
                    if (waitingQueue.size() > peakQueueLength) {
                        peakQueueLength = waitingQueue.size();
                        peakHour = currentTime;
                    }
                }
                
                // If customer might abandon (reneging)
                if (enableReneging) {
                    Event abandonment;
                    abandonment.type = Event::ABANDONMENT;
                    abandonment.time = currentTime + customers[customerId].maximumWaitTime;
                    abandonment.customerId = customerId;
                    eventQueue.push(abandonment);
                }
            }
            
            // Schedule next arrival if we haven't reached the customer limit
            if (customers.size() < numCustomers) {
                Event nextArrival;
                nextArrival.type = Event::ARRIVAL;
                nextArrival.time = currentTime + generateExponential(lambda, gen);
                nextArrival.customerId = 0;  // Will be assigned at arrival
                eventQueue.push(nextArrival);
            }
            
        } else if (event.type == Event::DEPARTURE) {
            // Process departure
            int customerId = event.customerId;
            int serverId = event.serverId;
            
            // Mark server as available
            serverBusy[serverId] = false;
            
            // Customer completes service
            Customer& customer = customers[customerId];
            if (!customer.abandoned) {
                completedCustomers++;
                totalRevenue += customer.revenue;
            }
            
            // Check if there are waiting customers
            if (!waitingQueue.empty()) {
                // Get next customer from queue
                int nextCustomerId = waitingQueue.front();
                waitingQueue.erase(waitingQueue.begin());
                
                // If this customer hasn't abandoned
                if (!customers[nextCustomerId].abandoned) {
                    // Assign to server
                    serverBusy[serverId] = true;
                    customers[nextCustomerId].serviceStartTime = currentTime;
                    double serviceTime = generateExponential(mu, gen);
                    customers[nextCustomerId].departureTime = currentTime + serviceTime;
                    
                    // Schedule departure
                    Event departure;
                    departure.type = Event::DEPARTURE;
                    departure.time = customers[nextCustomerId].departureTime;
                    departure.customerId = nextCustomerId;
                    departure.serverId = serverId;
                    eventQueue.push(departure);
                }
                // If customer abandoned, try next one in queue
                else {
                    // We already counted this customer as abandoned
                    // Try to find a non-abandoned customer
                    bool foundCustomer = false;
                    while (!waitingQueue.empty() && !foundCustomer) {
                        nextCustomerId = waitingQueue.front();
                        waitingQueue.erase(waitingQueue.begin());
                        
                        if (!customers[nextCustomerId].abandoned) {
                            // Found a non-abandoned customer
                            serverBusy[serverId] = true;
                            customers[nextCustomerId].serviceStartTime = currentTime;
                            double serviceTime = generateExponential(mu, gen);
                            customers[nextCustomerId].departureTime = currentTime + serviceTime;
                            
                            // Schedule departure
                            Event departure;
                            departure.type = Event::DEPARTURE;
                            departure.time = customers[nextCustomerId].departureTime;
                            departure.customerId = nextCustomerId;
                            departure.serverId = serverId;
                            eventQueue.push(departure);
                            
                            foundCustomer = true;
                        }
                    }
                }
            }
            
        } else if (event.type == Event::ABANDONMENT) {
            // Process abandonment
            int customerId = event.customerId;
            Customer& customer = customers[customerId];
            
            // Only process if customer is still waiting (not yet in service)
            if (!customer.abandoned && customer.serviceStartTime == 0) {
                customer.abandoned = true;
                abandonedCustomers++;
                lostRevenue += customer.revenue;
                
                // Note: We don't remove from waitingQueue here to avoid complexity
                // Instead, we'll check abandonment status when dequeuing
            }
        }
    }
    
    // Calculate performance metrics
    double totalWaitingTime = 0.0;
    double totalSystemTime = 0.0;
    int actualCompleted = 0;
    
    for (const auto& customer : customers) {
        if (!customer.abandoned && customer.departureTime > 0) {
            totalWaitingTime += (customer.serviceStartTime - customer.arrivalTime);
            totalSystemTime += (customer.departureTime - customer.arrivalTime);
            actualCompleted++;
        }
    }
    
    // Find peak hour
    auto peakIt = std::max_element(hourlyQueueCounts.begin(), hourlyQueueCounts.end(),
                                   [](const auto& p1, const auto& p2) { return p1.second < p2.second; });
    
    if (peakIt != hourlyQueueCounts.end()) {
        peakHour = peakIt->first;
        peakQueueLength = peakIt->second;
    }
    
    // Prepare statistics
    ServerStats stats;
    stats.avgWaitingTime = actualCompleted > 0 ? totalWaitingTime / actualCompleted : 0;
    stats.avgSystemTime = actualCompleted > 0 ? totalSystemTime / actualCompleted : 0;
    stats.utilization = totalBusyTime / (currentTime * numServers);
    stats.idleTime = 1.0 - stats.utilization;
    stats.maxQueueLength = maxQueueLength;
    stats.emptyQueueProbability = static_cast<double>(emptyQueueCount) / (completedCustomers + abandonedCustomers + emptyQueueCount);
    stats.totalRevenue = totalRevenue;
    stats.peakHour = peakHour;
    stats.peakQueueLength = peakQueueLength;
    stats.abandonedCustomers = abandonedCustomers;
    stats.lostRevenue = lostRevenue;
    
    // Print results if requested
    if (printDetails) {
        outFile << "======== COFFEE SHOP M/M/" << numServers << (enableReneging ? " QUEUE WITH RENEGING" : " QUEUE") << " SIMULATION RESULTS ========" << std::endl;
        outFile << "Simulation Parameters:" << std::endl;
        outFile << "  - Arrival Rate (λ): " << lambda << " customers/hour" << std::endl;
        outFile << "  - Service Rate (μ): " << mu << " customers/hour per server" << std::endl;
        outFile << "  - Number of Servers: " << numServers << std::endl;
        outFile << "  - Number of Customers: " << numCustomers << std::endl;
        if (enableReneging) {
            outFile << "  - Customer Patience Scale: " << patienceScale << " hour(s)" << std::endl;
        }
        outFile << "  - Total Simulation Time: " << currentTime << " hours" << std::endl << std::endl;
        
        outFile << "Performance Metrics:" << std::endl;
        outFile << std::fixed << std::setprecision(4);
        outFile << "  1. Average Waiting Time in Queue: " << stats.avgWaitingTime << " hours (" 
                << stats.avgWaitingTime * 60 << " minutes)" << std::endl;
        
        outFile << "  2. Average Time in System: " << stats.avgSystemTime << " hours (" 
                << stats.avgSystemTime * 60 << " minutes)" << std::endl;
        
        outFile << "  3. Utilization Factor: " << stats.utilization * 100 << "%" << std::endl;
        
        outFile << "  4. Idle Time of Baristas: " << stats.idleTime * 100 << "%" << std::endl;
        
        outFile << "  5. Maximum Queue Length: " << stats.maxQueueLength << " customers" << std::endl;
        
        outFile << "  6. Probability of an Empty Queue: " << stats.emptyQueueProbability * 100 << "%" << std::endl;
        
        outFile << "  7. Peak Hour Analysis:" << std::endl;
        outFile << "     - Peak Hour: Hour " << static_cast<int>(stats.peakHour) << std::endl;
        outFile << "     - Maximum Queue Length at Peak: " << stats.peakQueueLength << " customers" << std::endl;
        
        if (enableReneging) {
            outFile << "  8. Customer Abandonment:" << std::endl;
            outFile << "     - Total Abandoned Customers: " << stats.abandonedCustomers << std::endl;
            outFile << "     - Abandonment Rate: " << (static_cast<double>(stats.abandonedCustomers) / numCustomers) * 100 << "%" << std::endl;
            outFile << "     - Lost Revenue Due to Abandonment: $" << stats.lostRevenue << std::endl;
        }
        
        outFile << "  " << (enableReneging ? "9" : "8") << ". Revenue Analysis:" << std::endl;
        outFile << "     - Total Revenue: $" << stats.totalRevenue << std::endl;
        if (enableReneging) {
            outFile << "     - Potential Revenue (if no abandonment): $" << stats.totalRevenue + stats.lostRevenue << std::endl;
        }
        outFile << std::endl;
    }
    
    return stats;
}

int main() {
    // Parameters
    const double LAMBDA = 10.0;  // arrival rate (customers per hour)
    const double MU = 15.0;      // service rate (customers per hour)
    const int NUM_CUSTOMERS = 500;
    const double PATIENCE_SCALE = 0.25;   // average patience of 15 minutes (0.25 hour)

    // Random number generation
    std::random_device rd;
    std::mt19937 gen(rd());

    std::cout << "Starting simulation..." << std::endl;

    // Simulation without reneging
    std::ofstream outFileNoReneging("coffee_shop_simulation_results.txt");
    if (!outFileNoReneging) {
        std::cerr << "Error opening output file!" << std::endl;
        return 1;
    }
    
    outFileNoReneging << "==============================================" << std::endl;
    outFileNoReneging << "   COFFEE SHOP QUEUE SIMULATION - NO RENEGING" << std::endl;
    outFileNoReneging << "==============================================" << std::endl << std::endl;

    std::cout << "Running baseline simulation without reneging..." << std::endl;
    ServerStats baselineStatsNoReneging = simulateQueue(LAMBDA, MU, 1, NUM_CUSTOMERS, false, 0, true, outFileNoReneging, gen);
    
    std::cout << "Running service rate experiments without reneging..." << std::endl;
    ServerStats fasterStatsNoReneging = simulateQueue(LAMBDA, 20.0, 1, NUM_CUSTOMERS, false, 0, true, outFileNoReneging, gen);
    ServerStats slowerStatsNoReneging = simulateQueue(LAMBDA, 12.0, 1, NUM_CUSTOMERS, false, 0, true, outFileNoReneging, gen);
    
    std::cout << "Running multiple server simulation without reneging..." << std::endl;
    ServerStats multiServerStatsNoReneging = simulateQueue(LAMBDA, MU, 2, NUM_CUSTOMERS, false, 0, true, outFileNoReneging, gen);
    
    // Comparative Analysis
    outFileNoReneging << "==============================================" << std::endl;
    outFileNoReneging << "               COMPARATIVE ANALYSIS          " << std::endl;
    outFileNoReneging << "==============================================" << std::endl << std::endl;
    
    outFileNoReneging << std::fixed << std::setprecision(2);
    outFileNoReneging << "Average Waiting Time Comparison:" << std::endl;
    outFileNoReneging << "  - Baseline (M/M/1, μ = 15): " << baselineStatsNoReneging.avgWaitingTime * 60 << " minutes" << std::endl;
    outFileNoReneging << "  - Faster Service (M/M/1, μ = 20): " << fasterStatsNoReneging.avgWaitingTime * 60 << " minutes" << std::endl;
    outFileNoReneging << "  - Slower Service (M/M/1, μ = 12): " << slowerStatsNoReneging.avgWaitingTime * 60 << " minutes" << std::endl;
    outFileNoReneging << "  - Two Baristas (M/M/2, μ = 15 per server): " << multiServerStatsNoReneging.avgWaitingTime * 60 << " minutes" << std::endl;
    outFileNoReneging << std::endl;
    
    outFileNoReneging << "Utilization Factor Comparison:" << std::endl;
    outFileNoReneging << "  - Baseline (M/M/1, μ = 15): " << baselineStatsNoReneging.utilization * 100 << "%" << std::endl;
    outFileNoReneging << "  - Faster Service (M/M/1, μ = 20): " << fasterStatsNoReneging.utilization * 100 << "%" << std::endl;
    outFileNoReneging << "  - Slower Service (M/M/1, μ = 12): " << slowerStatsNoReneging.utilization * 100 << "%" << std::endl;
    outFileNoReneging << "  - Two Baristas (M/M/2, μ = 15 per server): " << multiServerStatsNoReneging.utilization * 100 << "%" << std::endl;
    outFileNoReneging << std::endl;
    
    outFileNoReneging.close();
    
    // Simulation with reneging
    std::ofstream outFileReneging("coffee_shop_simulation_with_reneging_results.txt");
    if (!outFileReneging) {
        std::cerr << "Error opening output file!" << std::endl;
        return 1;
    }
    
    outFileReneging << "==============================================" << std::endl;
    outFileReneging << "   COFFEE SHOP QUEUE SIMULATION WITH RENEGING" << std::endl;
    outFileReneging << "==============================================" << std::endl << std::endl;

    std::cout << "Running baseline simulation with reneging..." << std::endl;
    ServerStats baselineStatsReneging = simulateQueue(LAMBDA, MU, 1, NUM_CUSTOMERS, true, PATIENCE_SCALE, true, outFileReneging, gen);
    
    std::cout << "Running service rate experiments with reneging..." << std::endl;
    ServerStats fasterStatsReneging = simulateQueue(LAMBDA, 20.0, 1, NUM_CUSTOMERS, true, PATIENCE_SCALE, true, outFileReneging, gen);
    ServerStats slowerStatsReneging = simulateQueue(LAMBDA, 12.0, 1, NUM_CUSTOMERS, true, PATIENCE_SCALE, true, outFileReneging, gen);
    
    std::cout << "Running multiple server simulation with reneging..." << std::endl;
    ServerStats multiServerStatsReneging = simulateQueue(LAMBDA, MU, 2, NUM_CUSTOMERS, true, PATIENCE_SCALE, true, outFileReneging, gen);
    
    // Comparative Analysis
    outFileReneging << "==============================================" << std::endl;
    outFileReneging << "               COMPARATIVE ANALYSIS          " << std::endl;
    outFileReneging << "==============================================" << std::endl << std::endl;
    
    outFileReneging << std::fixed << std::setprecision(2);
    outFileReneging << "Average Waiting Time Comparison:" << std::endl;
    outFileReneging << "  - Baseline (M/M/1, μ = 15): " << baselineStatsReneging.avgWaitingTime * 60 << " minutes" << std::endl;
    outFileReneging << "  - Faster Service (M/M/1, μ = 20): " << fasterStatsReneging.avgWaitingTime * 60 << " minutes" << std::endl;
    outFileReneging << "  - Slower Service (M/M/1, μ = 12): " << slowerStatsReneging.avgWaitingTime * 60 << " minutes" << std::endl;
    outFileReneging << "  - Two Baristas (M/M/2, μ = 15 per server): " << multiServerStatsReneging.avgWaitingTime * 60 << " minutes" << std::endl;
    outFileReneging << std::endl;
    
    outFileReneging << "Abandonment Rate Comparison:" << std::endl;
    outFileReneging << "  - Baseline (M/M/1, μ = 15): " << (static_cast<double>(baselineStatsReneging.abandonedCustomers) / NUM_CUSTOMERS) * 100 << "%" << std::endl;
    outFileReneging << "  - Faster Service (M/M/1, μ = 20): " << (static_cast<double>(fasterStatsReneging.abandonedCustomers) / NUM_CUSTOMERS) * 100 << "%" << std::endl;
    outFileReneging << "  - Slower Service (M/M/1, μ = 12): " << (static_cast<double>(slowerStatsReneging.abandonedCustomers) / NUM_CUSTOMERS) * 100 << "%" << std::endl;
    outFileReneging << "  - Two Baristas (M/M/2, μ = 15 per server): " << (static_cast<double>(multiServerStatsReneging.abandonedCustomers) / NUM_CUSTOMERS) * 100 << "%" << std::endl;
    outFileReneging << std::endl;
    
    outFileReneging << "Lost Revenue Comparison:" << std::endl;
    outFileReneging << "  - Baseline (M/M/1, μ = 15): $" << baselineStatsReneging.lostRevenue << std::endl;
    outFileReneging << "  - Faster Service (M/M/1, μ = 20): $" << fasterStatsReneging.lostRevenue << std::endl;
    outFileReneging << "  - Slower Service (M/M/1, μ = 12): $" << slowerStatsReneging.lostRevenue << std::endl;
    outFileReneging << "  - Two Baristas (M/M/2, μ = 15 per server): $" << multiServerStatsReneging.lostRevenue << std::endl;
    outFileReneging << std::endl;
    
    outFileReneging << "Revenue Impact Analysis:" << std::endl;
    outFileReneging << "  - Additional revenue from faster service: $" << baselineStatsReneging.lostRevenue - fasterStatsReneging.lostRevenue << std::endl;
    outFileReneging << "  - Revenue loss from slower service: $" << slowerStatsReneging.lostRevenue - baselineStatsReneging.lostRevenue << std::endl;
    outFileReneging << "  - Revenue benefit from adding second barista: $" << baselineStatsReneging.lostRevenue - multiServerStatsReneging.lostRevenue << std::endl;
    outFileReneging << std::endl;
    
    outFileReneging.close();
    
    std::cout << "Simulation completed. Results saved to:" << std::endl;
    std::cout << "- coffee_shop_simulation_results.txt (without reneging)" << std::endl;
    std::cout << "- coffee_shop_simulation_with_reneging_results.txt (with reneging)" << std::endl;
    
    return 0;
}