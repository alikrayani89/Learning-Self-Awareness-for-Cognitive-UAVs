%% GNG function
function net = GrowingNeuralGasNetwork(inputData, params, plotFlag, GSVdimensions)

% alpha for derivative
alphaDist = 0.85;
alphaParam = alphaDist/(GSVdimensions/2);
betaDist = 1 - alphaDist;
betaParam = betaDist/(GSVdimensions/2);

wVec(1,1:(GSVdimensions/2)) = betaParam;
wVec(1,(GSVdimensions/2)+1:GSVdimensions) = alphaParam;

%% GNG parameters
N = params.N;                                                               % number of nodes in GNG
MaxIt = params.MaxIt;
L_growing = params.L_growing;
epsilon_b = params.epsilon_b;                                               % stepsize to update weight of winner node
epsilon_n = params.epsilon_n;                                               % stepsize to update weight of all direct neibours
alpha = params.alpha;
delta = params.delta;
T = params.T;
k = params.k;
L_decay = params.L_decay;

seedvector  = params.seedvector;

% Normalize input  data
minDataNorm = min(inputData);
dataNorm = inputData - repmat(minDataNorm,size(inputData,1),1);
maxDataNorm = max(dataNorm);
inputNorm = dataNorm./repmat(maxDataNorm,size(inputData,1),1);
inputNormOrd = inputNorm;

% Randomize input data
nData = size(inputNorm,1);                                                  % Size of input data (number of training samples)
nDim = size(inputNorm,2);                                                   % Dimension of input data
rng(seedvector)
inputNorm = inputNorm(randperm(nData), :);                                  % Random permutation of input data

CVmin = min(inputNorm);
CVmax = max(inputNorm);


% Initialization
Ni = 2;                                                                     % Initial 2 nodes for training the algorithm

wNorm = zeros(Ni, nDim);
for i = 1:Ni
    wNorm(i,:) = unifrnd(CVmin, CVmax);                                     % It returns an array of random numbers generated from the continuous uniform distributions with lower and upper endpoints specified by 'CVmin' and 'CVmax'.
end

E = zeros(Ni,1);
utility = ones(Ni,1);
C = zeros(Ni, Ni);
t = zeros(Ni, Ni);

% Loop

nx = 0;

for it = 1:MaxIt
    for c = 1:nData                                                         %   Number of data
        % Select Input
        nx = nx + 1;                                                        %   Counter of cycles inside the algorithm
        x = inputNorm(c,:);                                                 %   pick first input vector from permuted inputs
        
        % Competion and Ranking
        
        %         d = pdist2(x, wNorm,'euclidean');                         %   pdist- Pairwise distance between two sets of observations(Eucledian distance between input and 2 nodes initialised before)
        d = (wNorm - x).^2;
        d = repmat(wVec,size(d,1),1).*d;
        d = sqrt(sum(d,2));
        
        [~, SortOrder] = sort(d);                                           %   Organize distances between nodes and the first data point in an ascending order
        %         dbstop if error
        s1 = SortOrder(1);                                                  %   Closest node index to the first data point
        s2 = SortOrder(2);                                                  %   Second closest node index to the first data point
        
        % Aging of the closest neuron
        t(s1, :) = t(s1, :) + 1;                                            %   Increment the age of all edges emanating from s1
        t(:, s1) = t(:, s1) + 1;
        
        % Add Error
        dist0  = d(s1)^2;
        dist1  = d(s2)^2;
        E(s1) = E(s1) + dist0;
        % utility
        
        deltaUtility =  dist1 - dist0;                                      % Initial utility is zero in first case and dist is the error of first node
        utility(s1) =  utility(s1) + deltaUtility ;                         % Difference between error of two nodes
        
        
        % Adaptation
        wNorm(s1,:) = wNorm(s1,:) + epsilon_b*(x-wNorm(s1,:));              %   Move the nearest distance node and it's neibors towards input signal by fractions Eb and En resp.
        Ns1 = find(C(s1,:) == 1);                                           %   Take all the connections of the closest node to the data in question
        
        for j = Ns1
            wNorm(j,:) = wNorm(j,:) + epsilon_n*(x-wNorm(j,:));             %   Move the direct topological neibors of nearest distance node (S1) and it's neibors to wards input signal by fractions Eb and En resp.
        end
        
        % Create Link
        C(s1,s2) = 1;                                                       %  If s1 and s2 are connected by an edge , set the age of this edge to zero , it such edge doesn't exist create it
        C(s2,s1) = 1;
        t(s1,s2) = 0;                                                       %   Age of the edge
        t(s2,s1) = 0;
        % Remove Old Links
        C(t > T) = 0;                                                       %   remove edges with an age larger than Amax(a threshold value)
        nNeighbor = sum(C);                                                 %   Number of conecctions of each node
        AloneNodes = (nNeighbor==0);
        C(AloneNodes, :) = [];
        C(:, AloneNodes) = [];
        t(AloneNodes, :) = [];
        t(:, AloneNodes) = [];
        wNorm(AloneNodes, :) = [];
        E(AloneNodes) = [];
        utility(AloneNodes) = [];
        
        % Add New Nodes
        if mod(nx, L_growing) == 0 && size(wNorm,1) < N
            [~, q] = max(E);                                                %   Determine the unit q with the maximum accumulated error
            [~, f] = max(C(:,q).*E);                                        %   Maximum index related to the error related to a connected node
            
            r = size(wNorm,1) + 1;                                          %   Total number of nodes
            wNorm(r,:) = (wNorm(q,:) + wNorm(f,:))/2;                       %   Insert a new unit r halfway between q and it's neibor f with the largest error variable
            
            %   Remove old connections and introduce the presence of the
            %   new created node
            C(q,f) = 0;
            C(f,q) = 0;
            C(q,r) = 1;
            C(r,q) = 1;
            C(r,f) = 1;
            C(f,r) = 1;
            t(r,:) = 0;
            t(:, r) = 0;
            
            E(q) = alpha*E(q);                                              %   Decrease the error variable of q and f by multiplying them with a constand 'alpha'
            E(f) = alpha*E(f);
            %             E(r) = 0.5 *( E(q) + E(f) );                      %   Initialize the error of the new node equal to error of the winner node
            E(r) =   E(q) ;                                                 %   Initialize the error of the new node equal to error of the winner node
            
            utility(r) = 0.5 *( utility(q) + utility(f) );
        end
        
        % Eliminate Nodes
        if mod(nx, L_decay) == 0
            
            [max_E, ~] = max(E);                                            % Maximum accumelated error
            
            [min_utility,node_useless] = min(utility);                      % Node node_useless having minimum utility
            
            CONST = min_utility * k;                                        % Utility factor
            
            if (CONST < max_E)
                C(node_useless,:) = [];                                     % Remove the connection having smaller utility factor                                     % Remove the node having smaller utility factor
                C(:,node_useless) = [];
                wNorm(node_useless,:) = [];                                 % Remove the node having smaller utility factor
                utility(node_useless) = [];                                 % Remove the min utility value from the utility vector
                E(node_useless) = [];                                       % Remove error vector correspond to the node having min utility
                t(node_useless,:) = [];                                     % Remove agging vector correspond to the node having min utility
                t(:, node_useless) = [];
                
            end
            
        end  % end of decaying loop
        % Decrease Errors
        E = delta * E;                                                  % Decrease error variables by multiplying them with a constant delta
        utility = delta*utility;                                       % Decrease the utility by alpha_utility constant
    end     % end of data samples loop
    
    %% plot clustering results
    if plotFlag
        GlobalIteration = it;
        figure(1)
        
        inputNormObj = inputNorm(:,1:GSVdimensions);
        
        wNormObj = wNorm(:, 1:GSVdimensions);
        
        PlotResults(inputNormObj, wNormObj, C);
        pause(0.01);
    end
end  % end of iteration loop

%% Export Results data samples in nodes
datanodesNorm = cell(1,size(wNorm,1));
datanodes = cell(1,size(wNorm,1));
dataColorNode = [];

for c = 1:nData                                                             %   Number of data
    x = inputNormOrd(c,:);
    d = (wNorm - x).^2;
    d = repmat(wVec,size(d,1),1).*d;
    d = sqrt(sum(d,2));
    [~, minNode] = min(d);                                                  %   Organize distances between nodes and the first data point in an ascending order
    dataColorNode = [dataColorNode; minNode];
    datanodesNorm{1,minNode} = [datanodesNorm{1,minNode}; x];               % normalize and ordered data
    
    x = inputData(c,:);      %not normalize and orderd data
    datanodes{1,minNode} = [datanodes{1,minNode}; x];
end

%% Denormalization of GNG generated centroids
w = wNorm.*repmat(maxDataNorm,size(wNorm,1),1);
w = w + repmat(minDataNorm,size(wNorm,1),1);

%% >
indexAA = 0;
for indexA = 1:size(datanodes,2)
    if isempty(datanodes{1,indexA})
        % dont add it
    else
        indexAA = indexAA + 1;
        datanodes2{1,indexAA} = datanodes{1,indexA};
        datanodesNorm2{1,indexAA} = datanodesNorm{1,indexA};
    end
end
datanodes = datanodes2;
datanodesNorm = datanodesNorm2;
%% <


%%  Calculation of mean and covariances of generated nodes
nodesMean = [];
nodesMeanNorm = [];
for i = 1:size(datanodes,2)
    %   Calculation of mean values
    if size(datanodesNorm{1,i}, 1) > 1
        nodesMean = [nodesMean; mean(datanodes{1,i})];
        muX=nodesMean;
        assignin('base','muX',nodesMean)
        nodesMeanNorm = [nodesMeanNorm; mean(datanodesNorm{1,i})];
    else
        nodesMean = [nodesMean; (datanodes{1,i})]; 
        assignin('base','muX',nodesMean)
        nodesMeanNorm = [nodesMeanNorm; (datanodesNorm{1,i})]; 
    end 
    
    %   Calculation of covariance values
    nodesCov{1,i} = cov(datanodes{1,i});
    covX=nodesCov{1,i};
    assignin('base','covX',nodesCov{1,i})
    nodesCovNorm{1,i} = cov(datanodesNorm{1,i});
    
    if size(datanodesNorm{1,i}, 1) > 1
        nodesRadAccept_old(1,i) = sqrt(sum((3*std(datanodes{1,i}(:,1:(GSVdimensions/2)))).^2));        
        nodesRadAcceptNorm_old(1,i) =sqrt(sum((3*std(datanodesNorm{1,i}(:,1:(GSVdimensions/2)))).^2));
        
    else
        nodesRadAccept_old(1,i) = 0;     
        nodesRadAcceptNorm_old(1,i) =0;
    end
end

%% outputs of GNG
%   Data inside each node
net.datanodes  = datanodes;
net.datanodesNorm  = datanodesNorm ;

net.N = size(datanodes,2);  %   number of nodes
net.wNorm = wNorm;   % protocol values
net.w = w;   % protocol values (calculated centroids)
net.E = E;   % error value
net.C = C;   % connections
net.t = t;   % connection values
net.Parameters = params;
net.dataColorNode = dataColorNode;

%   Nodes' mean values and covariances
net.nodesMean = nodesMean;
net.nodesMeanNorm = nodesMeanNorm;
net.nodesCov= nodesCov;
net.nodesCovNorm = nodesCovNorm;


%   Nodes' radius of acceptances
net.nodesRadAccept_old = nodesRadAccept_old;
net.nodesRadAcceptNorm_old = nodesRadAcceptNorm_old;

net.data = inputNormOrd;

% Nodes' Mahalanobis Distance
% % net.distanceMaha = distanceMaha;

end