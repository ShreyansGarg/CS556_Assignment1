#include <iostream>
#include <queue>
#include <vector>
#include <random>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <algorithm>
#include <numeric>

struct Customer {
    double arrivalTime;
    double serviceStartTime;
    double serviceEndTime;
    double serviceTime;
};

struct Period {
    double startTime;
    double endTime;
    int numTellers;
    std::string name;
    
    // Statistics for each period
    int totalCustomers = 0;
    double totalWaitingTime = 0.0;
    double totalSystemTime = 0.0;
    double totalBusyTime = 0.0;
    double totalQueueLength = 0.0;
    int samplesForQueueLength = 0;
    int allTellersBusyCount = 0;
    int observations = 0;
    
    std::vector<double> queueLengthOverTime;
    std::vector<double> waitingTimeOverTime;
    std::vector<double> timePoints;
};

class BankSimulation {
private:
    double simulationEndTime;
    double arrivalRate;  // customers per hour
    double serviceTime;  // hours per customer
    
    std::vector<Period> periods;
    std::vector<Customer> completedCustomers;
    std::vector<std::pair<double, double>> queueLengthHistory;
    std::vector<std::pair<double, double>> waitTimeHistory;

    // Output file stream
    std::ofstream outputFile;

    // Random number generation
    std::mt19937 rng;
    std::exponential_distribution<double> arrivalDist;
    std::exponential_distribution<double> serviceDist;

public:
    BankSimulation(double simEndTime, double arrRate, double servTime, const std::string& outputFileName)
        : simulationEndTime(simEndTime), arrivalRate(arrRate), serviceTime(servTime) {
        
        // Open output file
        outputFile.open(outputFileName);
        if (!outputFile.is_open()) {
            std::cerr << "Error: Could not open output file " << outputFileName << std::endl;
            exit(1);
        }
        
        // Initialize random number generators with correct parameters
        std::random_device rd;
        rng = std::mt19937(rd());
        
        // For exponential distribution, parameter is lambda (rate parameter)
        // For interarrival times: lambda = arrivalRate
        arrivalDist = std::exponential_distribution<double>(arrivalRate);
        
        // For service times: lambda = 1/serviceTime (reciprocal of mean service time)
        serviceDist = std::exponential_distribution<double>(1.0 / serviceTime);
        
        // Define the time periods (in hours from start of day)
        periods.push_back({0.0, 2.0, 2, " 9:00 AM - 11:00 AM"});       // 9:00 AM - 11:00 AM: 2 tellers
        periods.push_back({2.0, 5.0, 4, "11:00 AM - 2:00 PM"});       // 11:00 AM - 2:00 PM: 4 tellers
        periods.push_back({5.0, 8.0, 3, "2:00 PM - 5:00 PM"});        // 2:00 PM - 5:00 PM: 3 tellers
    }
    
    ~BankSimulation() {
        if (outputFile.is_open()) {
            outputFile.close();
        }
    }
    
    // Generate time until next arrival (in hours)
    double generateInterArrivalTime() {
        return arrivalDist(rng);
    }
    
    // Generate service time for a customer (in hours)
    double generateServiceTime() {
        return serviceDist(rng);
    }
    
    int getTellersForTime(double time) {
        for (const auto &period : periods) {
            if (time >= period.startTime && time < period.endTime) {
                return period.numTellers;
            }
        }
        return 0; // Outside business hours
    }
    
    Period& getPeriodForTime(double time) {
        for (auto &period : periods) {
            if (time >= period.startTime && time < period.endTime) {
                return period;
            }
        }
        // Default to the first period (should not happen)
        return periods[0];
    }
    
    void runSimulation() {
        std::queue<Customer> customerQueue;
        std::vector<double> tellerEndTimes; // When each teller will be free
        tellerEndTimes.resize(4, 0.0); // Maximum of 4 tellers
        
        double currentTime = 0.0;
        double nextArrivalTime = currentTime + generateInterArrivalTime();
        
        double lastEventTime = currentTime;
        
        // Record queue length at the beginning
        queueLengthHistory.push_back({currentTime, 0.0});
        
        int customersArrived = 0;
        int customersServed = 0;
        
        while (currentTime < simulationEndTime) {
            // Update queue length stats for the time between last event and current time
            if (currentTime > lastEventTime) {
                double timeElapsed = currentTime - lastEventTime;
                
                // Update period-specific queue length
                Period& period = getPeriodForTime(lastEventTime);
                period.totalQueueLength += customerQueue.size() * timeElapsed;
                period.samplesForQueueLength++;
                
                if (period.timePoints.empty() || currentTime > period.timePoints.back() + 0.05) {
                    period.queueLengthOverTime.push_back(customerQueue.size());
                    period.timePoints.push_back(currentTime);
                }
                
                // Count busy tellers
                int busyTellers = 0;
                int availableTellers = getTellersForTime(lastEventTime);
                for (int i = 0; i < availableTellers; ++i) {
                    if (tellerEndTimes[i] > lastEventTime) {
                        busyTellers++;
                        period.totalBusyTime += std::min(tellerEndTimes[i], currentTime) - lastEventTime;
                    }
                }
                
                // Check if all tellers are busy at this time
                if (busyTellers == availableTellers && availableTellers > 0) {
                    period.allTellersBusyCount++;
                }
                period.observations++;
                
                lastEventTime = currentTime;
            }
            
            // Determine the next event (arrival or service completion)
            double nextServiceCompletion = simulationEndTime;
            int nextTellerToFinish = -1;
            
            // Find the next teller to finish service
            for (size_t i = 0; i < tellerEndTimes.size(); ++i) {
                if (tellerEndTimes[i] > currentTime && tellerEndTimes[i] < nextServiceCompletion) {
                    nextServiceCompletion = tellerEndTimes[i];
                    nextTellerToFinish = i;
                }
            }
            
            // Process the next event
            if (nextArrivalTime < nextServiceCompletion && nextArrivalTime < simulationEndTime) {
                // Process arrival
                currentTime = nextArrivalTime;
                customersArrived++;
                
                Customer newCustomer;
                newCustomer.arrivalTime = currentTime;
                newCustomer.serviceTime = generateServiceTime();
                
                // Check if any teller is available
                int numTellers = getTellersForTime(currentTime);
                int availableTeller = -1;
                
                for (int i = 0; i < numTellers; ++i) {
                    if (tellerEndTimes[i] <= currentTime) {
                        availableTeller = i;
                        break;
                    }
                }
                
                if (availableTeller != -1) {
                    // Teller is available, start service immediately
                    newCustomer.serviceStartTime = currentTime;
                    newCustomer.serviceEndTime = currentTime + newCustomer.serviceTime;
                    tellerEndTimes[availableTeller] = newCustomer.serviceEndTime;
                    
                    // Update statistics for this customer
                    double waitingTime = 0.0; // No waiting
                    Period& period = getPeriodForTime(currentTime);
                    period.totalWaitingTime += waitingTime;
                    period.totalSystemTime += waitingTime + newCustomer.serviceTime;
                    period.totalCustomers++;
                    period.waitingTimeOverTime.push_back(waitingTime);
                    
                    completedCustomers.push_back(newCustomer);
                    customersServed++;
                    
                    // Record waiting time
                    waitTimeHistory.push_back({currentTime, waitingTime});
                } else {
                    // All tellers busy, join the queue
                    customerQueue.push(newCustomer);
                }
                
                // Schedule next arrival
                nextArrivalTime = currentTime + generateInterArrivalTime();
                
                // Record queue length
                queueLengthHistory.push_back({currentTime, static_cast<double>(customerQueue.size())});
                
            } else if (nextServiceCompletion < simulationEndTime) {
                // Process service completion
                currentTime = nextServiceCompletion;
                tellerEndTimes[nextTellerToFinish] = 0.0; // Mark teller as free
                
                // If there are customers in the queue, serve the next one
                if (!customerQueue.empty()) {
                    Customer nextCustomer = customerQueue.front();
                    customerQueue.pop();
                    
                    nextCustomer.serviceStartTime = currentTime;
                    nextCustomer.serviceEndTime = currentTime + nextCustomer.serviceTime;
                    tellerEndTimes[nextTellerToFinish] = nextCustomer.serviceEndTime;
                    
                    double waitingTime = nextCustomer.serviceStartTime - nextCustomer.arrivalTime;
                    
                    // Update statistics for this customer
                    Period& period = getPeriodForTime(nextCustomer.arrivalTime);
                    period.totalWaitingTime += waitingTime;
                    period.totalSystemTime += waitingTime + nextCustomer.serviceTime;
                    period.totalCustomers++;
                    period.waitingTimeOverTime.push_back(waitingTime);
                    
                    completedCustomers.push_back(nextCustomer);
                    customersServed++;
                    
                    // Record waiting time
                    waitTimeHistory.push_back({currentTime, waitingTime});
                }
                
                // Record queue length
                queueLengthHistory.push_back({currentTime, static_cast<double>(customerQueue.size())});
            } else {
                // End of simulation
                currentTime = simulationEndTime;
            }
        }
        
        outputFile << "Customers arrived: " << customersArrived << std::endl;
        outputFile << "Customers served: " << customersServed << std::endl;
    }
    
    void printResults() {
        outputFile << "Bank Simulation Results" << std::endl;
        outputFile << "=======================" << std::endl;
        
        // Calculate overall statistics
        double totalWaitingTime = 0.0;
        double totalSystemTime = 0.0;
        int totalCustomers = 0;
        double totalBusyTime = 0.0;
        int allTellersBusyCount = 0;
        int totalObservations = 0;
        
        // Print period-by-period statistics
        for (const auto& period : periods) {
            outputFile << "\nPeriod: " << period.name << " (" << period.numTellers << " tellers)" << std::endl;
            outputFile << "-------------------------------------------------" << std::endl;
            
            outputFile << "Total customers in period: " << period.totalCustomers << std::endl;
            
            double avgWaitingTime = period.totalCustomers > 0 ? period.totalWaitingTime / period.totalCustomers : 0;
            double avgSystemTime = period.totalCustomers > 0 ? period.totalSystemTime / period.totalCustomers : 0;
            double avgQueueLength = period.totalQueueLength / (period.endTime - period.startTime);
            double utilizationRate = period.totalBusyTime / (period.numTellers * (period.endTime - period.startTime));
            double allTellersBusyProb = period.observations > 0 ? static_cast<double>(period.allTellersBusyCount) / period.observations : 0;
            
            outputFile << "a. Average waiting time: " << std::fixed << std::setprecision(2) << avgWaitingTime * 60 << " minutes" << std::endl;
            outputFile << "b. Average time in system: " << std::fixed << std::setprecision(2) << avgSystemTime * 60 << " minutes" << std::endl;
            outputFile << "c. Teller utilization rate: " << std::fixed << std::setprecision(2) << utilizationRate * 100 << "%" << std::endl;
            outputFile << "d. Average queue length: " << std::fixed << std::setprecision(2) << avgQueueLength << " customers" << std::endl;
            outputFile << "e. Probability of all tellers busy: " << std::fixed << std::setprecision(2) << allTellersBusyProb * 100 << "%" << std::endl;
            
            // Accumulate overall statistics
            totalWaitingTime += period.totalWaitingTime;
            totalSystemTime += period.totalSystemTime;
            totalCustomers += period.totalCustomers;
            totalBusyTime += period.totalBusyTime;
            allTellersBusyCount += period.allTellersBusyCount;
            totalObservations += period.observations;
        }
        
        // Print overall statistics
        outputFile << "\nOverall Statistics (9:00 AM - 5:00 PM)" << std::endl;
        outputFile << "-------------------------------------------------" << std::endl;
        outputFile << "Total customers served: " << totalCustomers << std::endl;
        
        double avgWaitingTime = totalCustomers > 0 ? totalWaitingTime / totalCustomers : 0;
        double avgSystemTime = totalCustomers > 0 ? totalSystemTime / totalCustomers : 0;
        
        // Calculate average queue length across all recorded points
        double totalQueueLength = 0.0;
        for (const auto& period : periods) {
            totalQueueLength += period.totalQueueLength;
        }
        double avgQueueLength = totalQueueLength / simulationEndTime;
        
        // Calculate overall utilization
        double overallUtilization = 0.0;
        double totalTellerHours = 0.0;
        for (const auto& period : periods) {
            totalTellerHours += period.numTellers * (period.endTime - period.startTime);
        }
        overallUtilization = totalBusyTime / totalTellerHours;
        
        double allTellersBusyProb = totalObservations > 0 ? static_cast<double>(allTellersBusyCount) / totalObservations : 0;
        
        outputFile << "a. Average waiting time: " << std::fixed << std::setprecision(2) << avgWaitingTime * 60 << " minutes" << std::endl;
        outputFile << "b. Average time in system: " << std::fixed << std::setprecision(2) << avgSystemTime * 60 << " minutes" << std::endl;
        outputFile << "c. Overall teller utilization rate: " << std::fixed << std::setprecision(2) << overallUtilization * 100 << "%" << std::endl;
        outputFile << "d. Average queue length: " << std::fixed << std::setprecision(2) << avgQueueLength << " customers" << std::endl;
        outputFile << "e. Probability of all tellers busy: " << std::fixed << std::setprecision(2) << allTellersBusyProb * 100 << "%" << std::endl;
        
        // Output data for plotting
        outputPlotData("queue_length.csv", queueLengthHistory);
        outputPlotData("waiting_time.csv", waitTimeHistory);
        
        // Analysis and recommendations
        outputFile << "\nAnalysis:" << std::endl;
        outputFile << "-------------------------------------------------" << std::endl;
        
        // Calculate theoretical utilization for each period
        outputFile << "Theoretical utilization rates:" << std::endl;
        for (const auto& period : periods) {
            // λ/(c*μ) where λ is arrival rate, μ is service rate, c is number of tellers
            double lambda = arrivalRate;                   // customers per hour
            double mu = 1.0 / serviceTime;                 // customers per hour
            double rho = lambda / (period.numTellers * mu);
            outputFile << "  " << period.name << ": " << std::fixed << std::setprecision(2) 
                      << rho * 100 << "%" << std::endl;
        }
        
        outputFile << "\nFindings:" << std::endl;
        bool hasInsufficientStaffing = false;
        for (size_t i = 0; i < periods.size(); ++i) {
            double lambda = arrivalRate;
            double mu = 1.0 / serviceTime;
            double rho = lambda / (periods[i].numTellers * mu);
            
            if (rho >= 1.0) {
                outputFile << "- The " << periods[i].name << " period has insufficient tellers, as the theoretical utilization exceeds 100%." << std::endl;
                hasInsufficientStaffing = true;
            }
        }
        
        if (!hasInsufficientStaffing) {
            outputFile << "- All periods have sufficient tellers based on theoretical utilization rates." << std::endl;
        }
        
        // Compare periods
        int worstPeriodIndex = 0;
        double maxWaitTime = 0;
        for (size_t i = 0; i < periods.size(); ++i) {
            double avgWaitTime = periods[i].totalCustomers > 0 ? periods[i].totalWaitingTime / periods[i].totalCustomers : 0;
            if (avgWaitTime > maxWaitTime) {
                maxWaitTime = avgWaitTime;
                worstPeriodIndex = i;
            }
        }
        
        outputFile << "- The period with the longest average waiting time is " << periods[worstPeriodIndex].name << "." << std::endl;
        
        // Recommendations
        outputFile << "\nRecommendations:" << std::endl;
        
        // For each period, calculate the minimum number of tellers needed to maintain stability
        for (size_t i = 0; i < periods.size(); ++i) {
            double lambda = arrivalRate;
            double mu = 1.0 / serviceTime;
            int minTellers = static_cast<int>(std::ceil(lambda / mu));
            if (periods[i].numTellers < minTellers) {
                outputFile << "- Increase the number of tellers during " << periods[i].name 
                          << " from " << periods[i].numTellers << " to at least " << minTellers << "." << std::endl;
            }
        }
        
        // Optimize teller allocation
        int currentTotal = 0;
        for (const auto& period : periods) {
            currentTotal += period.numTellers * (period.endTime - period.startTime);
        }
        
        outputFile << "- Proposed optimal teller allocation (maintaining total of " << currentTotal << " teller-hours):" << std::endl;
        
        // Simple reallocation based on observed waiting times
        std::vector<double> waitingTimeRatios;
        double totalWaitTimeRatio = 0.0;
        
        for (const auto& period : periods) {
            double avgWaitTime = period.totalCustomers > 0 ? period.totalWaitingTime / period.totalCustomers : 0;
            double duration = period.endTime - period.startTime;
            double waitTimeRatio = (avgWaitTime + 0.01) * duration; // Weighted by period duration
            waitingTimeRatios.push_back(waitTimeRatio);
            totalWaitTimeRatio += waitTimeRatio;
        }
        
        // Normalize ratios and calculate new allocation
        std::vector<int> newAllocation;
        double minTellerFraction = 0.0; // To store the min tellers needed across all periods
        
        for (size_t i = 0; i < periods.size(); ++i) {
            double lambda = arrivalRate;
            double mu = 1.0 / serviceTime;
            double minTellerNeeded = lambda / mu; // Minimum theoretical tellers needed
            minTellerFraction += minTellerNeeded * (periods[i].endTime - periods[i].startTime);
            
            double normalized = waitingTimeRatios[i] / totalWaitTimeRatio;
            double allocation = normalized * currentTotal / (periods[i].endTime - periods[i].startTime);
            int suggestedTellers = std::max(static_cast<int>(std::ceil(minTellerNeeded)), 
                                            static_cast<int>(std::round(allocation)));
            newAllocation.push_back(suggestedTellers);
        }
        
        // Adjust to match the current total teller-hours if needed
        int newTotal = 0;
        for (size_t i = 0; i < newAllocation.size(); ++i) {
            newTotal += newAllocation[i] * (periods[i].endTime - periods[i].startTime);
        }
        
        // If we're over budget, reduce tellers where possible
        while (newTotal > currentTotal) {
            int indexToReduce = -1;
            double minImpact = std::numeric_limits<double>::max();
            
            for (size_t i = 0; i < periods.size(); ++i) {
                double lambda = arrivalRate;
                double mu = 1.0 / serviceTime;
                double minRequired = lambda / mu;
                
                if (newAllocation[i] > std::ceil(minRequired)) {
                    double avgWaitTime = periods[i].totalCustomers > 0 ? 
                                        periods[i].totalWaitingTime / periods[i].totalCustomers : 0;
                    double duration = periods[i].endTime - periods[i].startTime;
                    
                    // Metric: lower wait time means less impact when reducing tellers
                    if (avgWaitTime < minImpact) {
                        minImpact = avgWaitTime;
                        indexToReduce = i;
                    }
                }
            }
            
            if (indexToReduce != -1) {
                newAllocation[indexToReduce]--;
                newTotal -= (periods[indexToReduce].endTime - periods[indexToReduce].startTime);
            } else {
                // Can't reduce further without going below minimum requirements
                break;
            }
        }
        
        // If we're under budget, add tellers where most needed
        while (newTotal < currentTotal) {
            int indexToIncrease = -1;
            double maxWaitTime = -1;
            
            for (size_t i = 0; i < periods.size(); ++i) {
                double avgWaitTime = periods[i].totalCustomers > 0 ? 
                                    periods[i].totalWaitingTime / periods[i].totalCustomers : 0;
                                    
                if (avgWaitTime > maxWaitTime) {
                    maxWaitTime = avgWaitTime;
                    indexToIncrease = i;
                }
            }
            
            if (indexToIncrease != -1) {
                newAllocation[indexToIncrease]++;
                newTotal += (periods[indexToIncrease].endTime - periods[indexToIncrease].startTime);
            } else {
                // No periods with waiting times (shouldn't happen)
                break;
            }
        }
        
        // Output the new allocation
        for (size_t i = 0; i < periods.size(); ++i) {
            outputFile << "  " << periods[i].name << ": " << periods[i].numTellers << " -> " << newAllocation[i] << " tellers" << std::endl;
        }
        
        // Additional recommendations
        outputFile << "\nAdditional recommendations:" << std::endl;
        outputFile << "- Consider implementing an appointment system for non-urgent transactions to distribute customer arrivals more evenly." << std::endl;
        outputFile << "- Evaluate cross-training employees to flexibly adjust staffing during unexpected peak times." << std::endl;
        outputFile << "- Monitor customer arrival patterns regularly to refine staffing models as patterns change." << std::endl;
    }
    
    void outputPlotData(const std::string& filename, const std::vector<std::pair<double, double>>& data) {
        std::ofstream outFile(filename);
        outFile << "Time,Value" << std::endl;
        
        for (const auto& point : data) {
            outFile << point.first << "," << point.second << std::endl;
        }
        
        outFile.close();
    }
};

int main() {
    // Simulation parameters - fixed with correct units
    double simulationEndTime = 8.0;        // 8-hour workday
    double arrivalRate = 40.0;             // 40 customers per hour
    double serviceTime = 4.0 / 60.0;       // 4 minutes = 4/60 hours per customer
   // std::ofstream outFile("output.txt"); // Create or open the file
    // Create and run the simulation
    BankSimulation simulation(simulationEndTime, arrivalRate, serviceTime,"report.txt");
    simulation.runSimulation();
    simulation.printResults();
    std::cout<<"the results are in report.txt";
    
    return 0;
}