% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                                           %
%    CRNtoMAS                                                               %
%                                                                           %
%                                                                           %
% OUTPUT: Returns the system of ordinary differential equations for a mass  %
%    action system. The output variables 'N', 'reactant_complexes',         %
%    'x_dot', 'ode', and 'model' allow the user to view the following,      %
%    respectively:                                                          %
%       - Stochiometric matrix                                              %
%       - Matrix of reaction complexes                                      %
%       - Left-hand side of the equations                                   %
%       - Right-hand side of the equations                                  %
%       - Complete network with all the species listed in the 'species'     %
%            field of the structure 'model'                                 %
%                                                                           %
% INPUT: model: a structure, representing the CRN (see README.txt for       %
%    details on how to fill out the structure)                              %
%                                                                           %
% Notes:                                                                    %
%    1. It is assumed that the network has mass action kinetics.            %
%    2. Ideas for some parts of the code was motivated by Soranzo and       %
%          Altafini.                                                        %
%                                                                           %
% Reference: Soranzo N, Altafini C (2009) ERNEST: a toolbox for chemical    %
%    reaction network theory. Bioinform 25(21):2853–2854.                   %
%    https://doi.org/10.1093/bioinformatics/btp513                          %
%                                                                           %
% Created: 24 June 2023                                                     %
% Last Modified: 28 June 2023                                               %
%                                                                           %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



function [N, reactant_complexes, x_dot, ode, model] = CRNtoMAS(model)

%
% Step 1: Create a list of all species indicated in the reactions
%

% Initialize list of species
model.species = { };

% Get all species from reactants
for i = 1:numel(model.reaction)
    for j = 1:numel(model.reaction(i).reactant)
        model.species{end+1} = model.reaction(i).reactant(j).species;
    end
end

% Get species from products
for i = 1:numel(model.reaction)
    for j = 1:numel(model.reaction(i).product)
        model.species{end+1} = model.reaction(i).product(j).species;
    end
end

% Get only unique species
model.species = unique(model.species);

% Use lowercase letters for concentration of corresponding species
concentration = lower(model.species);

% Convert the concentration letters to symbolic form
species_concentration = sym(concentration, 'real');



%
% Step 2: Form stoichiometric matrix N
%

% Count the number of species
m = numel(model.species);

% Initialize the matrix of reactant complexes
reactant_complexes = [ ];

% Initialize the matrix of product complexes
product_complexes = [ ];

% Initialize the stoichiometric matrix
N = [ ];

% For each reaction in the model
for i = 1:numel(model.reaction)
  
    % Initialize the vector for the reaction's reactant complex
    reactant_complexes(:, end+1) = zeros(m, 1);
    
    % Fill it out with the stoichiometric coefficients of the species in the reactant complex
    for j = 1:numel(model.reaction(i).reactant)
        reactant_complexes(find(strcmp(model.reaction(i).reactant(j).species, model.species), 1), end) = model.reaction(i).reactant(j).stoichiometry;
    end
    
    % Initialize the vector for the reaction's product complex
    product_complexes(:, end+1) = zeros(m, 1);
    
    % Fill it out with the stoichiometric coefficients of the species in the product complex
    for j = 1:numel(model.reaction(i).product)
        product_complexes(find(strcmp(model.reaction(i).product(j).species, model.species), 1), end) = model.reaction(i).product(j).stoichiometry;
    end
    
    % Create a vector for the stoichiometric matrix: Difference between the two previous vectors
    N(:, end+1) = product_complexes(:, end) - reactant_complexes(:, end);
    
    % If the reaction is reversible
    if model.reaction(i).reversible
      
        % Insert a new vector for the reactant complex: make it same as the product complex
        reactant_complexes(:, end+1) = product_complexes(:, end);
        
        % Insert a new vector for the product complex: make it the same as the reactant complex
        product_complexes(:, end+1) = reactant_complexes(:, end-1);
        
        % Insert a new vector in the stoichiometric matrix: make it the additive inverse of the vector formed earlier
        N(:, end+1) = -N(:, end);
    end
end

% Count the total number of reactions
r = size(N, 2);



%
% Step 3: Form the reaction equations
%

% Create vector of symbolic rate constants
k = sym(strcat('k', string(1:r)), 'real');

% Create matrix of species
species_matrix = repmat(species_concentration', 1, r);

% Place corresponding species in reactant_complexes matrix
reactant_complexes_species = species_matrix.^reactant_complexes;

% Get the actual reactant complexes
reactant_complexes_vector = prod(reactant_complexes_species', 2);

% Multiply the rate constant with its corresponding complex
v = diag(k)*reactant_complexes_vector;



%
% Step 4: Form the mass action system
%

ode = N*v;



%
% Step 5: Display the mass action system
%

% Create vector of derivatives
x_dot = sym(strcat(concentration, '_dot'), 'real');

% Create cell array of derivatives, equal signs, and reactions
mass_action_system = horzcat(string(x_dot'), repmat('=',m, 1), string(ode));

% Get the variable letter used as species name
species_variable = regexprep(concentration,'[^a-zA-Z\s]','');

% Remove these from the species concentration variables
species_number_as_string = erase(concentration, species_variable);

% Display result depending on format of species used
try

    % For X1, X2, etc. format, convert strings to numbers
    species_number_as_number = cellfun(@str2num, species_number_as_string);

    % Rearrange the species
    [~, species_order] = sort(species_number_as_number);

    % Print the result
    fprintf('The mass action system for %s is: \n\n', model.id)
    for i = 1:m
        disp(strjoin(mass_action_system(species_order(i),:)))
    end
catch

    % Otherwise, just print the result without arranging
    fprintf('The mass action system for %s is: \n\n', model.id)
    for i = 1:m
        disp(strjoin(mass_action_system(i,:)))
    end
end

end