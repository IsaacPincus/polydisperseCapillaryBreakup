function [zet, varargout]=zetacal(Q,d,phieq,nk,hsk,himodel)
% Outputs:
% zeta - friction coefficient rescaled by Zimm friction coefficient (friction of
%        isolated, isotropic coils at equilibrium)
% G     - (optional) correction factor to Rouse-like drag of a chain of blobs
% zb    - (optional) inner-most blob's friction coefficient
% Nb    - (optional) number of inner-most blob's
%
% Inputs:
% Q    - rms end-to-end distance scaled by equilibirum value Q0
% d    - transverse diameter of partially stretched coil scaled by Q0; 
%        irrelevant parameter at present; set to 1
% phieq - equilibrium coil volume fraction = c/c*
% nk    - Number of Kuhn segments
% hsk - hydrodynamic interaction parameter of a Kuhn segment
% himodel - 0 or 1 - Constant drag calculated as the Rouse drag of a chain
%                    of Nk Kuhn segments with a segmental friction corresponding to the given
%                    hstar parameter
%           2 - Conformation-dependent drag without any concentration dependence
%           3 - Conformation-dependent drag including concentration dependence
%
% Calls:
% ubcal - friction coefficient of an isolated, isotropic random-walk scaled by 3 pi etas Q0

if hsk == 0 
    hsk = eps; 
end

e1 = exp(1);
L = nk^0.5;

switch himodel
    case 0 % constant friction coefficient is the same as the freely draining Rouse drag; no HI
        zet = 2*pi^0.5*hsk*nk^0.5/ubcal(hsk,nk); 
        G = 1;  % irrelevant in this model
        zb = 1; % irrelevant in this model
        Nb = 1; % irrelevant in this model
    
    case 1 % constant friction coefficient equal to fully stretched chain with concentration-dependent Batchelor correction        
           % Same as Case 3 below applied for fully stretched chains        
           xip = 1/nk^0.5; % t-blob size for fully-stretched chain
           Np = nk;        % no. of t-blobs for fully-tretched chains
           A = nk;         % Aspect ratio of t-blob rod for fully-tretched chains
           Q = nk^0.5;     % fully stretched chain
           % Single t-blob drag
           zb = ubcal(hsk,nk/Np)*xip/ubcal(hsk,nk);  % scaled by Zimm drag
           
           % Rouse-like drag of chain of t-blobs
           zeta_Rouse = zb*Np;
           
           % Correction due to inter-blob hydrodynamic interactions
           J = 1 + (0.5*(A-1) + (2/e1)*(A-1)*(A-1))/A;
           if phieq > 0 && phieq < nk^0.5
               phipincus = phieq*(Q*xip^2); % volume fraction of P-blob rods of length Q and dia xip
               chibyxip = (log(1/phipincus)/phipincus)^0.5; % hydrodynamic screening length for P-blob rods rods - ref. Mackaplow & Shaqfeh, JFM (1996), uol. 329, pp. 155-186
               B = chibyxip + 1;
               G1      = 1+ log( (J*(B-1) + A)/( B + A  )   ); % inspired by Batchelor's interpolation
               %                G      = zeta_Rouse*exp(-(Q-1)/2) + (1-exp(-(Q-1)/2))*log( (J*(X-1) + A)/( X + A  )   ); % inspired by Batchelor's interpolation
               G2 = 1 + log( (  (B-1) + 1)/( B + 1/A)   );   % correcting factor to achieve corrn = 1 in certain limiting cases
               
               G = G1/G2;
           elseif phieq == 0
               G1 = 1 + log( J);
               G2 = 1;%log( J  );
               G = G1/G2;
           else
               G = 1;
           end
           
           zet = zeta_Rouse/G;
           Nb = nk;
        
    case 2 % Only stretch dependence - no concentration dependence
        
        % Equilibrium coil drag
        % Nb = 1; Nk per blob = Nk
        zeta_isolcoil = ubcal(hsk,nk); % scaled by 3 pi etas Q0
        
        % Dividing by zeta_isolcoil  gives zeta normalized by Zimm friction coefficient
        
        % Partially stretched chain drag
        xip = max(1/Q,1/nk^0.5); % t-blob size
        Np = min(Q^2,nk);  % No. of t-blobs
        A = Q/xip; % Aspect ratio of t-blob pole
        
        % Single t-blob drag
        Ub = ubcal(hsk,nk/Np); % U_b = zeta_b/(3 pi etas xi_b) - for single t-blob
        zb = Ub*xip;  % t-blob drag rescaled by 3 pi etas Q0
        zb = zb/zeta_isolcoil;   % t-blob drag rescaled by Zimm friction coefficient
        Nb = Np; %Inner blobs
        
        % Rouse-like drag of chain of t-blobs
        zeta_Rouse = zb*Np;
        
        % Correction due to inter-blob hydrodynamic interactions
        J = 1 + (0.5*(A-1) + (2/e1)*(A-1)*(A-1))/A;
        G1 = 1+log(J);
        G2 = 1;
        G = G1/G2;
        
        zet = zeta_Rouse/G;
        
        if zet < 1
            zet = 1;
        end
        
    case 3 % Stretch & concentration dependence
 
        % Equilibrium coil drag
        % Nb = 1; Nk per blob = Nk
        zeta_isolcoil = ubcal(hsk,nk); % scaled by 6 pi etas Q0
        % Dividing by zeta_isolcoil  gives zeta normalized by Zimm friction coefficient
        
        % Coiled-state friction coefficient
        if phieq > L
            % concentrated (unentangled) regime
            Ub = ubcal(hsk,1); % U_b = zeta_K/(3 pi etas b_k) - for single Kuhn segment
            zb = Ub*1/sqrt(nk)/zeta_isolcoil; % Kuhn-step drag rescaled by by Zimm friction coefficient
            zeta_coil = zb*nk; % Rouse-like drag of chain of Kuhn segments; limiting value for Kramer's chain
            Nc = nk;
        elseif phieq > 1
            % Single c-blob drag
            Nc = phieq^2; % No. of c-blobs
            xic = 1/phieq; % c-blob size
            Ub = ubcal(hsk,nk/Nc); % U_b = zeta_b/(3 pi etas xi_b) - for single c-blob
            zb = Ub*xic;  % c-blob drag rescaled by 6 pi etas Q0
            zb = zb/zeta_isolcoil;   % c-blob drag rescaled by Zimm friction coefficient
            % Rouse-like drag of chain of c-blobs
            zeta_coil = zb*Nc; 
        else
            zeta_coil = 1; % Zimm friction coefficient
            Nc = 1;
            zb = 1;
        end
        
        % Partially stretched chain drag
        if Q < nk^0.5
            xip = 1/Q; % t-blob size
            Np = Q^2;  % No. of t-blobs
            A = Q/xip; % Aspect ratio of t-blob rod
            
            if phieq < Q % weak-screening of partially stretched chains
                Nb = Np;
                % Single t-blob drag
                Ub = ubcal(hsk,nk/Np); % U_b = zeta_b/(3 pi etas xi_b) - for single P-blob
                zb = Ub*xip;  % t-blob drag rescaled by 6 pi etas Q0
                zb = zb/zeta_isolcoil;   % t-blob drag rescaled by Zimm friction coefficient
                
                % Rouse-like drag of chain of t-blobs
                zeta_Rouse = zb*Np;
                
                % Correction due to inter-blob hydrodynamic interactions
                J = 1 + (0.5*(A-1) + (2/e1)*(A-1)*(A-1))/A;
                if phieq > 0
                    phipincus = phieq*(Q*xip^2); % volume fraction of P-blob rods of length Q and dia xip
                    chibyxip = (log(1/phipincus)/phipincus)^0.5; % hydrodynamic screening length for P-blob rods rods - ref. Mackaplow & Shaqfeh, JFM (1996), uol. 329, pp. 155-186
                    B = chibyxip + 1;
                    G1      = 1+ log( (J*(B-1) + A)/( B + A  )   ); % inspired by Batchelor's interpolation
                    %                G      = zeta_Rouse*exp(-(Q-1)/2) + (1-exp(-(Q-1)/2))*log( (J*(X-1) + A)/( X + A  )   ); % inspired by Batchelor's interpolation
                    G2 = 1 + log( (  (B-1) + 1)/( B + 1/A)   );   % correcting factor to achieve corrn = 1 in certain limiting cases
                    
                    G = G1/G2;
                else
                    G1 = 1 + log( J);
                    G2 = 1;%log( J  );
                    G = G1/G2;
                end
                
            else % strong-screening regime
                % No inter-blob hydrodynamic interactions i.e. corrn = 1
                %
                Nb = Nc;
                zeta_Rouse = zeta_coil;
                G = 1;
            end
            
        else % Gaussian chains can stretch beyond contour length; when that happens drag must be the drag for fully-stretched FE chains
             % Calculation same as for himodel = 1
             
            xip = 1/nk^0.5; % t-blob size for fully-stretched chain
            Np = nk;        % no. of t-blobs for fully-tretched chains
            Nb = Np;
            A = nk;         % Aspect ratio of t-blob rod for fully-tretched chains
            Q = nk^0.5;     % fully stretched chain
            % Single t-blob drag
            zb = ubcal(hsk,nk/Np)*xip/ubcal(hsk,nk);  % scaled by Zimm drag
            
            % Rouse-like drag of chain of t-blobs
            zeta_Rouse = zb*Np;
            
            % Correction due to inter-blob hydrodynamic interactions
            J = 1 + (0.5*(A-1) + (2/e1)*(A-1)*(A-1))/A;
            if phieq > 0 && phieq < nk^0.5
                phipincus = phieq*(Q*xip^2); % volume fraction of P-blob rods of length Q and dia xip
                chibyxip = (log(1/phipincus)/phipincus)^0.5; % hydrodynamic screening length for P-blob rods rods - ref. Mackaplow & Shaqfeh, JFM (1996), uol. 329, pp. 155-186
                B = chibyxip + 1;
                G1      = 1+ log( (J*(B-1) + A)/( B + A  )   ); % inspired by Batchelor's interpolation
                %                G      = zeta_Rouse*exp(-(Q-1)/2) + (1-exp(-(Q-1)/2))*log( (J*(X-1) + A)/( X + A  )   ); % inspired by Batchelor's interpolation
                G2 = 1 + log( (  (B-1) + 1)/( B + 1/A)   );   % correcting factor to achieve corrn = 1 in certain limiting cases
                
                G = G1/G2;
            elseif phieq == 0
                G1 = 1 + log( J);
                G2 = 1;%log( J  );
                G = G1/G2;
            else % phieq > nk^0.5
                G = 1;
            end
        end
            
        zet = zeta_Rouse/G;
        
%         if zet < zeta_coil
%             zet = zeta_coil;
%         end
        
end

zet = abs(zet);
varargout{1} = G;
varargout{2} = zb;
varargout{3} = Nb;
end


 






    