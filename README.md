# Power Control and Link Sceduling

The various parameters considered are obtained from the gains.csv and params.csv file. The gains.csv file contains the gain between every user pair, while the params.csv file consists of parameters such as :  
    n - Number of wireless users  
    theta - SINR Threshold for decoding  
    C - Weighting Parameter (Used to compute the objective function)  
    N - Noise Power  
In addition to selecting transmission power, I am also varying link activations making the system a TDMA one. This allows for a better performance in some cases.  
Two greedy approaches are implemented to find a good power allocation scheme along with a good TDMA scheme for link activation.
An objective function is used to gauge the performance of the algorithm. The higher the value of the objective function, the better its performance
