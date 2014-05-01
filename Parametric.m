% Solves parametric max-flow instances with 2 parameters.
%
%   E(x) = \sum_i (a_i+b_i\lambda + c_i\mu)x_i + \sum_(i,j\in N) x_i(1-x_j)*d_ij
%
% Note: the code only accept integer weights.

classdef Parametric < handle
	properties
		% Search space
		lambda_limit = [-1e3 1e3];
		mu_limit = [-1e3 1e3];
		
		iterations = 5000;
		
		solutions;
		
		plot_3d = false;
	end
	
	properties (SetAccess = protected)
		unary;
		slope1;
		slope2;
		
		pairwise;
	end
	
	methods (Static)
		
		% Simple and slow 2D implementation
		function E = create_connectivty(rows, cols, nbh)
			
			E = zeros(rows*cols*size(nbh,1),2);
			max_index = [rows cols];
			
			counter = 0;
			
			for ii = 1:rows;
				for jj = 1:cols;
					
					curr_index_linear = Parametric.linear_index([ii jj], max_index);
					
					for itr = 1:size(nbh,1)
						new_index = nbh(itr,:) + [ii jj];
						
						if (Parametric.valid_index(new_index, max_index))
							
							counter = counter+1;
							E(counter,:) = [curr_index_linear Parametric.linear_index(new_index, max_index)];
							
						end
					end
				end
			end
			
			% Purge
			for itr = 1:size(E,1);
				if (E(itr,1) == 0 && E(itr,2) == 0)
					startPurge = itr;
					break;
				end
			end
			
			E(startPurge:end,:) = [];
		end
		
		function valid = valid_index(index, max_index)
			valid = true;
			
			if any(index < 1)
				valid = false;
			end
			
			if any(index > max_index)
				valid = false;
			end
		end
		
		function lin_ind = linear_index(index, max_index)
			lin_ind = index(1) + (index(2)-1)*max_index(1);
		end
		
		
		function draw_solution_diagram(solutions, plot_3d)
			
			if nargin < 2
				plot_3d = false;
			end
			
			clf; hold on;
			
			for itr = 1:numel(solutions)
				if ( ~isempty(solutions))
					
					X = solutions(itr).polygon(1,:);
					Y = solutions(itr).polygon(2,:);
					Z = solutions(itr).polygon(3,:); %#ok<NASGU>
					
					drawnow()
					
					color = rand(1,3);
					
					if (plot_3d)
						fill3(X,Y,Z,color);
						view([-5 -7 3]);
					else
						fill(X,Y,color);
					end
				else
					disp 'Empty polygon';
				end
			end
			
			xlabel('lambda');
			ylabel('mu');
			title(sprintf(' %d solutions. \n Each polygon marks a region where a certain solution is optimal.', length(solutions)));
		end
		
	end
	
	methods
		
		function self = Parametric(unary, slope1, slope2, pairwise)
			% unary: 		A n x 1 matrix of all unary costs
			% slope1: 	A n x 1 matrix of the cost associated with lambda for each node.
			% slope2: 	A n x 1 matrix of the cost associated with mu for each node.
			% pairwise: A m x 4 of every pairwise cost formated as
			%            [node_1 node_2 cost_1 cost_2]
			%           where cost_1 is the cut cost from node_1 to node_2
			%           and vice versa.
			addpath([fileparts(mfilename('fullpath')) filesep 'include']);
			
			if (norm(unary - round(unary)) > 1e-10)
				error('Unary cost must be integers');
			end
			
			if (norm(slope1 - round(slope1)) > 1e-10)
				error('slope1 cost must be integers');
			end
			
			if (norm(slope2 - round(slope2)) > 1e-10)
				error('slope2 cost must be integers');
			end
			
			if (norm(pairwise(:,3:4) - round(pairwise(:,3:4))) > 1e-10)
				error('Pairwise costs must be integers');
			end
			
			% Check input
			if (size(unary) ~= size(slope1))
				error('Unary and slope1 must be same size');
			end
			
			if (size(slope1) ~= size(slope2))
				error('Slope1 and Slope2 must be same size');
			end
			
			if any(pairwise(:) < 0)
				error('Only positive weights');
			end
			
			if (size(pairwise,2) ~= 4)
				error('Pairwise should be a 4 columns matrix')
			end
			
			if max(max(pairwise(:,1:2))) > numel(unary)
				error('Invalid pairwise matrix');
			end
			
			if min(min(pairwise(:,1:2))) < 0
				error('Invalid pairwise matrix');
			end
			
			% Save input
			self.unary = round(unary);
			self.slope1 = round(slope1);
			self.slope2 = round(slope2);
			
			
			% Takes longer time to weed out duplicates then simply adding them again.
			pairwise(:,1:2) = pairwise(:,1:2)-1; % c++ is zero indexed.
			self.pairwise = pairwise;
		end
		
		
		
		% Output:
		%    A struct with k solutions. For each k.
		%
		%			solutions(k).polygon:			 The corners of the polygon defining where the
		%																 solutions is optimal [x y z];
		%			solutions(k).solution:		 The opitmal solutions.
		%			solutions(k).tanget_plane: The plane equation for this solutions given on affine form.
		function solutions = find_all_solutions(self)
			
			% compile
			cpp_file = 'parametric_mex.cpp';
			out_name = 'parametric_mex';
			
			extra_arguments = {};

			
			if ispc()
				% Set these paths (observe build paths and visual studio verion).
				% Note that the CGAL dll must be found by MATLAB. This can be done in two ways
				% 1) Append envoirment variable path to point to <cgal-build-folder>\bin 
				% 2) Copy it to the same folder as the mex file 
				boost_root = 'C:\dev\boost_1_55_0';
				cgal_root = 'C:\dev\CGAL-4.4';
				
				% Paths
				extra_arguments{end+1} = ['-I""' boost_root '""'];
				

				extra_arguments{end+1} = ['-I""' cgal_root '\include""'];
				extra_arguments{end+1} = ['-I""' cgal_root '\build\include""'];
				extra_arguments{end+1} = ['-I""' cgal_root '\auxiliary\gmp\include""'];
				
				extra_arguments{end+1} = ['-L""' boost_root '\lib64-msvc-12.0""'];
				extra_arguments{end+1} = ['-L""' cgal_root '\build\lib""'];
				extra_arguments{end+1} = ['-L""' cgal_root '\auxiliary\gmp\lib'];			
				extra_arguments{end+1} = ['-L""' cgal_root '\auxiliary\gmp\bin'];
				
				% Libaries (update here too)
				extra_arguments{end+1} = '-lCGAL-vc120-mt-4.4';
				extra_arguments{end+1} = '-llibboost_thread-vc120-mt-1_55';
				extra_arguments{end+1} = '-llibboost_system-vc120-mt-1_55';
				extra_arguments{end+1} = '-llibgmp-10';
				extra_arguments{end+1} = '-llibmpfr-4';
				
			else
				extra_arguments{end+1} = '-lCGAL';
				extra_arguments{end+1} = '-lboost_thread';
				extra_arguments{end+1} = '-lut';
				extra_arguments{end+1} = '-lmpfr';
				extra_arguments{end+1} = '-lgmp';
				
				extra_arguments{end+1} = '-I/usr/include';
				extra_arguments{end+1} = '-frounding-math';
			end	

			sources = {['maxflow-v3.02.src' filesep 'graph.cpp'], ...
				['utils' filesep 'cppmatrix.cpp'], ...
				['maxflow-v3.02.src' filesep 'graph.cpp'], ...
				[filesep 'maxflow-v3.02.src' filesep 'maxflow.cpp'], ...
				};
			
			compile(cpp_file, out_name, sources, extra_arguments);
			
			% Call solver
			interval = [self.lambda_limit self.mu_limit];
			solutions = parametric_mex(self.unary, self.slope1, self.slope2, self.pairwise, ...
				interval, self.iterations);
			
			self.solutions = solutions;
		end
		
		function display(self)
			details(self);
			
			if ~isempty(self.solutions)
				Parametric.plot_solution(self.solutions, self.plot_3d);
			end
		end
		
		
		function [labelling, energy] = evaluate(self, lambda,mu)
			% Evaluate a point (lambda,mu)
			
			cpp_file = 'single_point_mex.cpp';
			out_name = 'single_point_mex';
			extra_arguments{1} = ['-lut -lmpfr -lgmp -lboost_thread  -frounding-math'];
			
			sources = {['maxflow-v3.02.src' filesep 'graph.cpp'], ...
				['utils' filesep 'cppmatrix.cpp'], ...
				['maxflow-v3.02.src' filesep 'graph.cpp'], ...
				[filesep 'maxflow-v3.02.src' filesep 'maxflow.cpp'], ...
				};
			
			compile(cpp_file, out_name, sources, extra_arguments);
			
			[labelling,energy] = single_point_mex(self.unary, self.slope1, self.slope2, self.pairwise, [lambda mu]);
		end
	end
end