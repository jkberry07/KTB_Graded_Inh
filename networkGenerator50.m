%Generates a network of MC and GCs (15:1 GC:MC ratio) and saves data needed for simulation (note: distance_GC50.mat is calculated using granulecells50.mat)

%This code is mostly the same as that by Kersen et al., found at https://github.com/dkersen/olfactory-bulb
%(D. E. C. Kersen, G. Tavoni, and V. Balasubramanian. Connectivity and dynamics in the olfactory bulb. PLoS Comput Biol, 18(2):e1009856, 2022. ISSN 1553-7358.
%doi: 10.1371/journal.pcbi.1009856.)

%Calculated actual synapse location in addition to just the distance from MC in order to calculate distance 
%from "electrode" for LFP (original code on Github from Kersen 2022 didn't calculate the distance3D matrix) (Lines 76-78, 90, 113-114, 135-139)

EPLvolume = 1.5e+9; % volume of EPL in micrometers cubed (Richard et al 2010)
EPLthickness = 131; % measured thickness of EPL in micrometers
MCLthickness = 36;
IPLthickness = 27; 
EPLsurfarea = EPLvolume/EPLthickness; % calculates the surface assuming SA of inner and outer surfaces are approximately equal
totalGlomNum = 1800; % total number of glomeruli 
glomDensity = totalGlomNum/EPLsurfarea; % density of glomeruli

% maximum radius of the OB space, modeling the space as a thin cylinder
rmax = 425 %gives an area 50% the size of the original model

% number of glomeruli based on the area of the OB space
glomNum = round(glomDensity*(pi*rmax^2));

% number of MCs per glomerulus
glomitralNum = randi(11, 1, glomNum) + 14; % ~15-25 mitral cells per glomerulus 

% total number of mitral cells 
mitralNum = sum(glomitralNum); 

% generate the array of mitral cell objects
mitralArray(1:mitralNum) = mitral();

% matrix to hold the x-y locations of the glomeruli
glomXYarray = zeros(glomNum,3);

% index for assigning properties for MCs
currentMitralTotal = 0;
for i = 1:glomNum
    d = rmax + 1;
    
    % determine the xy-location of the glomeruli 
    r = rmax*sqrt(rand);
    theta = 2*pi*rand;
    glomerX = r*cos(theta);
    glomerY = r*sin(theta);
   
    % record the xy-location of the glomeruli
    glomXYarray(i,:) = [glomerX glomerY i];
    
    % for the MCs belonging to the given glomerulus, assign properties
    for j = 1:glomitralNum(i)
        mitralArray(currentMitralTotal+j) = mitralArray(currentMitralTotal+j).assignProperties(i, glomerX, glomerY, rmax, currentMitralTotal+j);
    end
    
    % update the MC index
    currentMitralTotal = currentMitralTotal + glomitralNum(i);
end

% set the target number of GCs
granPerMit = 15; 
granuleNum = granPerMit * mitralNum;

% create the array to store GC objects
granuleArray = [];

% create an empty matrix to record connections (a connnection between an
% MC and GC is recorded as 1 in the relevant entry in the network matrix,
% is 0 otherwise)
network = zeros(mitralNum, granuleNum);

% GC index 
granuleindex = 0;

% distance matrix for synapses during simulation (
distance = zeros(mitralNum, granuleNum);
% distance matrix for distance of synapse from electrode
distance_3D = zeros(mitralNum,granuleNum);
re = [0, 0,IPLthickness + MCLthickness + EPLthickness/2]; %position of electrode

while length(granuleArray) < granuleNum
    
    % generate a new GC and assign properties
    newGranule = granule();
    newGranule = newGranule.assignProperties(rmax);
    
    % generate an empty matrix representing potential connections between
    % the MCs of the network and the new GC
    tempNet = zeros(mitralNum,1);
    tempDistance = ones(mitralNum,1)*-1;
    tempDistance_3D = ones(mitralNum,1)*-1;


    % for each MC, determine the probability of synapse with the new GC
    mitralOrder = randperm(mitralNum); %Randomize order of MCs tested
    for j = 1:mitralNum
         
        % set the MC
        indx = mitralOrder(j);
        mitralCell = mitralArray(indx);
        
        % calculate the total pre-existing synapses 
        deadSpace = sum(network(j,:));
        
        % determine the location of the center of the GC cone at the height
        % of the MC
        [g_loc_x, g_loc_y] = newGranule.calculateLocation(mitralCell.z);
        
        % calculate the radius of the GC cone at the height of the MC
        g_radius = newGranule.calculateRadius(mitralCell.z);
        
        % calculate the distance between the GC cone and MC at the height
        % of the MC
        svec = [mitralCell.x, mitralCell.y] - [g_loc_x, g_loc_y];
        s = norm(svec);
        
        % if the MC and GC are in range of one another, calculate the
        % average number of synapses
        if s < mitralCell.radius + g_radius && mitralCell.z > newGranule.z0 && mitralCell.z < newGranule.zmax
            
            % the average number of synapses for the MC/GC pair
            lambda = synProb(newGranule, mitralCell, deadSpace, g_radius, s);
            
            % calculate the probability of synapse assuming a Poisson
            % distribution
            prob = 1 - exp(-lambda);
            
            % sample the given probability and assign synaptic distance
            if rand < prob
                tempNet(j,1) = 1;
                while tempDistance(j,1) == -1
                   r = rand*mitralCell.radius;
                   theta = rand * 2 * pi;
                   x_r = r*cos(theta);
                   y_r = r*sin(theta);
                   if norm([x_r,y_r] - svec) < g_radius 
                       tempDistance(j,1) = r;
                       %calculate distance_3D also
                       tempDistance_3D(j,1) = norm([x_r + mitralCell.x, y_r + mitralCell.y, mitralCell.z]-re);
                   end            
                end           
            end          
        end                 
    end                                          
    
    % calculate the total number of synapses the GC has made and
    % incorporate the new GC if it is greater than 0
    totalSynapses = sum(tempNet(:,1));
   
    if totalSynapses > 0
        % if the total number of synapses is greater than the number of
        % available spines, select a random subset of MCs equal to the
        % number of available spines to be incorporated into the network
        if totalSynapses > newGranule.availableSpines
           synmitrals = find(tempNet);
           synmitrals = synmitrals(randperm(length(synmitrals)));
           tempNet(:) = 0;
           tempNet(synmitrals(1:newGranule.availableSpines)) = 1;
           tempDistance(synmitrals(newGranule.availableSpines+1:totalSynapses)) = -1;
           tempDistance_3D(synmitrals(newGranule.availableSpines+1:totalSynapses)) = -1;
        end
        % update the GC index
        granuleindex = granuleindex + 1;
        
        % assign the new connections to the appropriate entry in the
        % network for the GC
        network(:, granuleindex) = tempNet;
        distance(:,granuleindex) = tempDistance;
        distance_3D(:,granuleindex) = tempDistance_3D;
        
        % update the GC array
        granuleArray = [granuleArray newGranule];
    end
    
    disp(length(granuleArray));
end

glomArray = zeros(1,mitralNum);
for i = 1:mitralNum
    glomArray(i) = mitralArray(i).glomerulus;
end

% save arrays and network matrix
save('mitralCells50.mat', 'mitralArray','-v7.3');
save('fullNetwork50.mat', 'network', '-v7.3');
save('granuleCells50.mat', 'granuleArray', '-v7.3');
save('glomeruli50.mat', 'glomXYarray','glomArray','-v7.3');
save('distance50.mat', 'distance','-v7.3');
save('distance_3D50.mat', 'distance_3D','-v7.3');
