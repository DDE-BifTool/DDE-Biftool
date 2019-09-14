function [newbranch, success] = br_bifdet(funcs, branch)
%% Compute stability, detect and correct bifurcations along branch
% function [newbranch, success] = br_bifdet(branch)
% Purpose:
%   This routine computes stability, detects possible bifurcations,
%   corrects bifurcation candidates to genuine bifurcation points,
%   and adds them to the branch as a "native" point of the branch
% INPUT:
%  funcs: defining functions
%	branch: continued branch (contains method, initial points and all new points)
% OUTPUT:
%	newbranch: extended branch
%  success: 1 if a point was detected, 0 otherwise
%
% $Id: br_bifdet_old.m 309 2018-10-28 19:02:42Z jansieber $
%
%%
last_sign = zeros(1,10);
last_sign(bif2num('genh')) = NaN;

if isempty(branch) % Use this to reset last_sign
   newbranch = branch;
   success = 1;
   return;
end

%% Process method parameters
detection_minimal_real_part = branch.method.bifurcation.minimal_real_part;
tolcorr = branch.method.bifurcation.correction_tolerance;
rtolfactor = branch.method.bifurcation.radial_tolerance_factor;
max_iter = branch.method.bifurcation.secant_iterations; % max number of secant iterations
conv_r = branch.method.bifurcation.secant_tolerance; % required accuracy for convergence
imagthresh=branch.method.bifurcation.imagthreshold; % threshold for treating a root as complex
isimag=@(x)x>imagthresh;
isreal=@(x)~isimag(x);
%% Set variables
success = 0;
kind = branch.point(1).kind;
newbranch = branch;
tempbranch = br_bifinit(branch); % Make sure all points get flag, nmfm, nvec (hopf)
method = branch.method.point;
stmethod = branch.method.stability;
ll = length(branch.point);
free_par = branch.parameter.free;

%% Check that we're on a supported type of branch
if ~(strcmp(kind,'stst') || strcmp(kind,'hopf'))
   error('BR_BIFDET: only detection on steady state and hopf branches is implemented.');
end

%% Ensure normal form, stability, null vectors are present in hopf branch
if strcmp(kind, 'hopf')
   if ~isfield(tempbranch.point(1).nmfm,'L1') || isempty(tempbranch.point(1).nmfm.L1) ||...
           ~isfield(tempbranch.point(1), 'nvec') || isempty(tempbranch.point(1).nvec)
      tempbranch.point(1) = nmfm_hopf(funcs,tempbranch.point(1));
   end
   if ~isfield(tempbranch.point(1),'stability') || isempty(tempbranch.point(1).stability)
      stability = p_stabil(funcs,tempbranch.point(1),branch.method.stability);
      tempbranch.point(1).stability = stability;
   end
end
%% Traverse branch
curind = 1;

while curind < length(tempbranch.point)
   
   curind = curind + 1;
   
   bifkind = '';
      
   new_point = tempbranch.point(curind);
   
   %% Make sure we have stability info
   if ~isfield(new_point,'stability') || isempty(new_point.stability)
      stability = p_stabil(funcs,new_point,branch.method.stability);
      if isempty(stability) || isempty(stability.l1)
         fprintf('BR_BIFDET: could not obtain roots of point %g.\n', curind);
         return;
      end
      new_point.stability = stability;
   else
      stability = new_point.stability;
   end
   roots = stability.l1;
   
   %% Make sure we have L1 and nullvectors in case of hopf branch 
   if strcmp(kind, 'hopf')
      if ~isfield(new_point.nmfm,'L1') || isempty(new_point.nmfm.L1) ||...
              ~isfield(new_point, 'nvec') || isempty(new_point.nvec)
         new_point = nmfm_hopf(funcs,new_point, tempbranch.point(curind-1));
      end
   end
   
   tempbranch.point(curind) = new_point;
   
   %% Select a subset of the roots to avoid overflow
   roots = roots(real(roots) >= detection_minimal_real_part);
   
   biffound = 0;
   continue_detection = 1;
   detection_status = zeros(1,10);
   
   %% Try to detect all bifurcations. 
   % If a bifurcation is detected, but fails to be corrected, continue to
   % detect the other bifurcations.
   while continue_detection
      
      biffound = 0;
      bifkind = 'NONE';
      
      %% Note: update this if more detections are added
      if strcmp(kind, 'stst')
         if sum(detection_status) >= 1 % We detect 1 bifurcation on stst branches (hopf)
            continue_detection = 0;
            break;
         end
      elseif strcmp(kind, 'hopf')
         if sum(detection_status) >= 3 % We detect 3 bifurcations on hopf branches
            continue_detection = 0;
            break;
         end
      end
      
      %% HOPF
      if detection_status(bif2num('hopf')) == 0 && strcmp(kind,'stst')
         new_sign = nmfm_hopfdet(roots);
         if new_sign ~= 0
            if abs(new_sign + last_sign(bif2num('hopf'))) < 1e-6
               % Hopf found!
               biffound = 1;
               bifkind = 'hopf';
            end
         end
         last_sign(bif2num('hopf')) = new_sign;
         detection_status(bif2num('hopf')) = 1;
      end
      
      %% FOLD
      if ~biffound && strcmp(kind,'stst')
         % Do fold detection
      end
      
      %% DOUBLE HOPF
      if ~biffound && detection_status(bif2num('hoho')) == 0 && strcmp(kind,'hopf')
         % Remove known roots, look for smallest real part of imaginary pair
         new_sign = nmfm_smrp(funcs, new_point,stmethod,...
             'remove_omega',true,'threshold',isimag);
         new_sign = new_sign/abs(new_sign);
         if abs(new_sign + last_sign(bif2num('hoho'))) < 1e-6
            % Double Hopf found!
            biffound = 1;
            bifkind = 'hoho';
         end
         last_sign(bif2num('hoho')) = new_sign;
         detection_status(bif2num('hoho')) = 1;
      end
      
      %% ZERO HOPF
      if ~biffound && detection_status(bif2num('zeho')) == 0 && strcmp(kind,'hopf')
         % Check whether the smallest real eigenvalue crosses 0
         new_sign = nmfm_smrp(funcs, new_point, stmethod,...
             'remove_omega',true,'threshold',isreal);
         new_sign = new_sign/abs(new_sign);
         if abs(new_sign + last_sign(bif2num('zeho'))) < 1e-6
            biffound = 1;
            %fprintf('BR_BIFDET: zero hopf detected at par(%g) = %g,
            %par(%g) = %g.\n', free_par(1),
            %new_point.parameter(free_par(1)), free_par(2),
            %new_point.parameter(free_par(2)));
            bifkind = 'zeho';
         end
         last_sign(bif2num('zeho')) = new_sign;
         detection_status(bif2num('zeho')) = 1;
      end
      
      %% GENERALIZED HOPF
      if ~biffound && detection_status(bif2num('genh')) == 0 && strcmp(kind,'hopf')
         % Do double hopf detection
         %new_point = nmfm_hopf(funcs,new_point);
         new_sign = new_point.nmfm.L1;
         new_sign = new_sign/abs(new_sign);
         if new_sign ~= last_sign(bif2num('genh')) && ~isnan(last_sign(bif2num('genh')))
            biffound = 1;
            bifkind = 'genh';
         end
         last_sign(bif2num('genh')) = new_sign;
         detection_status(bif2num('genh')) = 1;
      end
      
      %% NO MATCH
      if ~biffound
         break;
      end
      
      %% Construct point halfway between the new and old point
      halfway_point = p_axpy(1,new_point,tempbranch.point(curind-1));
      halfway_point = p_axpy(1/2,halfway_point,[]);
      
      if strcmp(bifkind, 'genh') || strcmp(bifkind, 'zeho') || strcmp(bifkind, 'hoho')
         % Keep new point as starting point for correction
         bifcand = new_point;
      else
         % Use halfway point as bifurcation candidate
         secant=p_axpy(-1,new_point,tempbranch.point(curind-1));
         secant=p_secant(secant,p_norm(new_point));
         bifcand = p_correc(funcs, halfway_point, free_par, secant, method);
         % Use this approximate stability only to convert to Hopf
         bifcand.stability = stability;
      end
      
      %% Change the new point to a bifurcation point
      good_conversion = 1;
      switch bifkind
         case 'hopf'
            if length(bifcand.stability.l1) < 2
               % Can not convert to Hopf point
               fprintf('BR_BIFDET: detected point has no suitable pair of eigenvalues.\n');
               good_conversion = 0;
            else
               new_bifpoint = p_tohopf(funcs,bifcand);
            end
         case 'fold'
            new_bifpoint = p_tofold(funcs,bifcand);
         case {'genh','hoho','zeho'}
            % Do nothing yet
         otherwise
            % This shouldn't ever happen because of earlier checks
            display(roots);
            error('BR_BIFDET: an unknown type of bifurcation was detected.\n');
      end
      
      if ~good_conversion
         fprintf('BR_BIFDET: Could not convert to %s point at par(%d) = %4.4f.\n',...
             bifkind, free_par(1), new_point.parameter(free_par(1)));
         continue;
      end
      
      %% Use secant methods for codim 2 bifurcations
      if strcmp(bifkind,'genh')
         %% Use secant method
         corrsuccess = 0;
         hopfpoint = bifcand;
         previous = tempbranch.point(curind-1);
         if hopfpoint.nmfm.L1 > 0
            pospoint = hopfpoint;
            negpoint = previous;
         else
            pospoint = previous;
            negpoint = hopfpoint;
         end
         prevpoint = previous;
         for i = 1:max_iter
            sumpoint = p_axpy(1, pospoint, negpoint);
            halfpoint = p_axpy(0.5, sumpoint, []);
            halfpoint = nmfm_hopf(funcs,halfpoint, prevpoint);
            %fprintf('P_CORREC: iteration %g: L1 = %g.\n', i, halfpoint.L1);
            if abs(halfpoint.nmfm.L1) < conv_r % done!
               corrsuccess = 1;
               break
            end
            if halfpoint.nmfm.L1 < 0
               negpoint = halfpoint;
            else
               pospoint = halfpoint;
            end
            prevpoint = halfpoint;
         end
         if corrsuccess
            new_bifpoint = p_togenh(halfpoint);
         end
      elseif strcmp(bifkind,'zeho')
         %% Use secant method
         corrsuccess = 0;
         hopfpoint = bifcand;
         hopfpoint.stability = p_stabil(funcs,bifcand,stmethod);
         prevpoint = tempbranch.point(curind-1);
         currp = nmfm_smrp(funcs, hopfpoint, stmethod,...
             'remove_omega',true,'threshold',isreal);
         if currp > 0
            pospoint = hopfpoint;
            negpoint = prevpoint;
         else
            pospoint = prevpoint;
            negpoint = hopfpoint;
         end
         for i = 1:max_iter
            sumpoint = p_axpy(1, pospoint, negpoint);
            halfpoint = p_axpy(0.5, sumpoint, []);
            [currp, stability] = nmfm_smrp(funcs, halfpoint, stmethod,...
                'remove_omega',false,'threshold',isreal);
            if abs(currp) < conv_r % Smallest real part zero
                halfpoint.stability=stability;
                corrsuccess = 1;
                break;
            end
            if currp < 0
               negpoint = halfpoint;
            else
               pospoint = halfpoint;
            end
         end
         if corrsuccess
            new_bifpoint = p_tozeho(halfpoint);
            last_sign(bif2num('genh'))=NaN;
         end
      elseif strcmp(bifkind,'hoho')
         %% Use secant method
         corrsuccess = 0;
         hopfpoint = bifcand;
         prevpoint = tempbranch.point(curind-1);
         cursign = nmfm_smrp(funcs, hopfpoint, stmethod,...
             'remove_omega',true,'threshold',isimag);
         if cursign > 0
            pospoint = hopfpoint;
            negpoint = prevpoint;
         else
            pospoint = prevpoint;
            negpoint = hopfpoint;
         end
         secant=p_secant(p_axpy(-1,negpoint,pospoint),p_norm(pospoint));
         for i = 1:max_iter
            sumpoint = p_axpy(1, pospoint, negpoint);
            halfpoint = p_axpy(0.5, sumpoint, []);
            halfpoint = p_correc(funcs, halfpoint,free_par,secant,method);
            [cursign, stability] = nmfm_smrp(funcs, halfpoint, stmethod,...
                'remove_omega',true,'threshold',isimag);
            if abs(cursign) < conv_r
               halfpoint.stability = stability;
               corrsuccess = 1;
               break;
            end
            if cursign < 0
               negpoint = halfpoint;
            else
               pospoint = halfpoint;
            end
         end
         if corrsuccess
            halfpoint = nmfm_hopf(funcs, halfpoint, prevpoint);
            new_bifpoint = p_tohoho(halfpoint);
         end
      else
          %% Correct codim 1 point
          % Change convergence criteria
          method.minimal_accuracy = tolcorr;
          [new_bifpoint, corrsuccess] = p_correc(funcs,new_bifpoint,free_par,[],method);
      end
      
      if ~corrsuccess
         fprintf('BR_BIFDET: Unable to correct as %s point near par(%d) = %4.4f.\n',...
             bifkind, free_par(1), new_point.parameter(free_par(1)));
         continue_detection = 1;
      else
         continue_detection = 0;
      end
      
   end % while loop for detection
   
   % Check whether loop actually produced a result
   if ~biffound
      continue;
   end
   
   %% From here on there is an actual bifurcation
   lenpar = length(free_par);
   if lenpar == 2
      fprintf('BR_BIFDET: %s point found at par(%d) = %.10f, par(%d) = %.10f.\n',...
          bifkind, free_par(1), new_bifpoint.parameter(free_par(1)),...
          free_par(2), new_bifpoint.parameter(free_par(2)));
   else
      fprintf('BR_BIFDET: %s point found at par(%d) = %.10f.\n',...
          bifkind, free_par(1), new_bifpoint.parameter(free_par(1)));
   end
   
   %% Convert to degenerate branch point
   switch kind
      case 'stst'
         new_degpoint = p_tostst(funcs,new_bifpoint);
      case 'hopf'
         if ~isfield(new_bifpoint, 'stability') || ~isfield(new_bifpoint.stability,'l1')...
                 || isempty(new_bifpoint.stability.l1)
            new_bifpoint.stability = p_stabil(funcs, new_bifpoint, branch.method.stability);
         end
         new_degpoint = p_tohopf(funcs,new_bifpoint);
      otherwise
         display(kind)
         error('BR_BIFDET: detection for this branch type is not supported.\n');
   end
   %% Region check
   % Check whether the point is in an acceptable region
   % The acceptable region is a cylinder with radius
   % rtolfactor*d(point(curind-1),point(curind))
   old_point = tempbranch.point(curind-1);
   oldtonew = p_axpy(-1, old_point, new_point); % new_point - old_point
   direction = p_axpy(1/p_norm(oldtonew),oldtonew,[]); % normalize
   oldtodeg = p_axpy(-1, old_point, new_degpoint); % new_degpoint - old_point
   projection = p_axpy(p_inprod(oldtodeg, direction),direction,old_point);
   oldtoproj = p_axpy(-1,old_point,projection);
   newtoproj = p_axpy(-1, oldtoproj, oldtodeg);
   % Hack to disregard v
   newtoproj.v = 0;
   radial_distance = p_norm(newtoproj);
   radiustol = rtolfactor*p_norm(oldtonew);
   
   line_distance = p_norm(p_axpy(-1,halfway_point, projection));
   linetol = p_norm(oldtonew)/2;
   
   if ~(line_distance <= linetol && radial_distance <= radiustol)
      fprintf('BR_BIFDET: the detected %s point does not fall within the branch.\n', bifkind);
      continue;
   end
   
   %% Normal form computation!
   switch bifkind
      case 'hopf'
         new_bifpoint = nmfm_hopf(funcs, new_bifpoint,  tempbranch.point(curind-1));
         fprintf('BR_BIFDET: L1 = %.10f, omega = %.10f, par(%d) = %.10f.\n',...
             new_bifpoint.nmfm.L1, new_bifpoint.omega, free_par(1),...
             new_bifpoint.parameter(free_par));
      case 'genh'
         new_bifpoint = nmfm_genh(funcs, new_bifpoint,  tempbranch.point(curind-1));
         fprintf('BR_BIFDET: L2 = %.10f, omega = %.10f, par(%d) = %.10f, par(%d) = %.10f.\n',...
             new_bifpoint.nmfm.L2, new_bifpoint.omega, free_par(1),...
             new_bifpoint.parameter(free_par(1)), free_par(2),...
             new_bifpoint.parameter(free_par(2)));
      case 'zeho'
         fprintf('BR_BIFDET: omega = %.10f, par(%d) = %.10f, par(%d) = %.10f.\n',...
             new_bifpoint.omega, free_par(1), new_bifpoint.parameter(free_par(1)),...
             free_par(2), new_bifpoint.parameter(free_par(2)));
         new_bifpoint = nmfm_zeho(funcs, new_bifpoint, tempbranch.point(curind-1));
      case 'hoho'
         fprintf('BR_BIFDET: omega1 = %.10f, omega2 = %.10f, par(%d) = %.10f, par(%d) = %.10f.\n',...
             new_bifpoint.omega1, new_bifpoint.omega2, free_par(1),...
             new_bifpoint.parameter(free_par(1)), free_par(2),...
             new_bifpoint.parameter(free_par(2)));
         new_bifpoint = nmfm_hoho(funcs, new_bifpoint, tempbranch.point(curind-1));
      otherwise
         % This shouldn't ever happen because of earlier checks
         display(bifkind);
         error('BR_BIFDET: an unknown type of bifurcation was detected.\n');
   end
   %% Create point structure fro mcomputated information
   % Set flag
   new_degpoint.flag = bifkind;
   
   % Add normal form information to branchpoint
   new_degpoint.nmfm = new_bifpoint.nmfm;
   new_degpoint.nvec = new_bifpoint.nvec;
   
   % Match presence of stability
   if ~isfield(tempbranch.point(curind-1),'stability')
      new_degpoint = rmfield(new_degpoint,'stability');
   else
      if ~isfield(new_degpoint,'stability')
         new_degpoint.stability = [];
      end
   end
   
   % Match presence of null vectors
   if ~isfield(tempbranch.point(curind-1),'nvec')
      new_degpoint = rmfield(new_degpoint,'nvec');
   else
      if ~isfield(new_degpoint,'nvec')
         new_degpoint.nvec = [];
      end
   end
   
   % Put in place between the two previous points
   tempbranch.point = [tempbranch.point(1:curind-1), new_degpoint, tempbranch.point(curind:ll)];
   
   % Increase current index once more because of insertion
   curind = curind + 1;
   
end

success = 1;
newbranch = tempbranch;
end
