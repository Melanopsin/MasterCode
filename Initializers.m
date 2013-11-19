close all
clear all
clc



% % function Melanopsin()
% % % MODEL OF MELANOPSIN ACTIVATION
% % 
% % % M0* + G.GDP <-- kG1/kG2 --> M0*.G.GDP
% % % M0*.G.GDP -- kG3 --> M0*.G + GDP
% % % M0*.G + GTP -- kG4 --> M0*.G.GTP
% % % M0*.G.GTP -- kG5 --> M0* + Ga.GTP + Gbg
% % % PLC + Ga.GTP -- kP --> PLC*.Ga.GTP
% % % PLC*.Ga.GTP -- kI1 --> PLC.Ga.GDP
% % % PLC.Ga.GDP -- kI2 --> PLC + Ga.GDP
% % % Ga.GDP + Gbg -- kI3 --> G.GDP
% % % PIP2 + PLC*.Ga.GTP -- kS --> SecM + PLC*.Ga.GTP
% % % SecM + Channel- -- kC --> SecM.Channel+
% % 
% % kG1 = 0;
% % kG2 = 0;
% % kG3 = 0;
% % kG4 = 0;
% % kG5 = 0;
% % kp = 0;
% % kI1 = 0;
% % kI2 = 0;
% % kI3 = 0;
% % kS = 0;
% % kC = 0;
% % 
% % %% SPECIES: X = [
% % %% X(1)           G.GDP 
% % %% X(2)           G.GTP 
% % %% X(3)           Ga.GTP
% % %% X(4)           Gbg
% % %% X(5)           PLC
% % %% X(6)           PLC*.Ga.GTP
% % %% X(7)           PLC.Ga.GDP
% % %% X(8)           Ga.GDP
% % %% X(9)           G.GDP
% % %% X(10)          PIP2
% % %% X(11)          SecM
% % %% X(12)          Channel-
% % %% X(13)          SecM.Channel+ ];
% % 
% % %% MELANOPSIN COMPLEXES: M = [
% % %% M(1)                        M0*
% % %% M(2)                        M0*.G.GDP
% % %% M(3)                        M0*.G
% % %% M(4)                        M0*.G.GTP
% % %% M(5)                        M0*.K
% % %% M(6)                        M1*
% % %% M(7)                        M1*.G.GDP
% % %% M(8)                        M1*.G
% % %% M(9)                        M1*.G.GTP
% % %% M(10)                       M1*.K
% % %% M(11)                       M1*.ArrB1
% % %% M(12)                       M1*.ArrB2
% % %% M(13)                       M2*
% % %% M(14)                       M2*.G.GDP
% % %% M(15)                       M2*.G
% % %% M(16)                       M2*.G.GTP
% % %% M(17)                       M2*.K
% % %% M(18)                       M2*.ArrB1
% % %% M(19)                       M2*.ArrB2
% % %% M(20)                       M3*
% % %% M(21)                       M3*.G.GDP
% % %% M(22)                       M3*.G
% % %% M(23)                       M3*.G.GTP
% % %% M(24)                       M3*.K
% % %% M(25)                       M3*.ArrB1
% % %% M(26)                       M3*.ArrB2
% % %% M(27)                       M4*
% % %% M(28)                       M4*.G.GDP
% % %% M(29)                       M4*.G
% % %% M(30)                       M4*.G.GTP
% % %% M(31)                       M4*.K
% % %% M(32)                       M4*.ArrB1
% % %% M(33)                       M4*.ArrB2
% % %% M(34)                       M5*
% % %% M(35)                       M5*.G.GDP
% % %% M(36)                       M5*.G
% % %% M(37)                       M5*.G.GTP
% % %% M(38)                       M5*.K
% % %% M(39)                       M5*.ArrB1
% % %% M(40)                       M5*.ArrB2
% % %% M(41)                       M6*
% % %% M(42)                       M6*.G.GDP
% % %% M(43)                       M6*.G
% % %% M(44)                       M6*.G.GTP
% % %% M(45)                       M6*.K
% % %% M(46)                       M6*.ArrB1
% % %% M(47)                       M6*.ArrB2 ];
% % 
% % 
% % 

dataset = 0;

% case 0 --- calcium imaging data


switch dataset
    
    case 0
        
        load('calciumimaging.mat')
        
    otherwise
        
end


% % 
% % t = 0; X = zeros(1,12); M = zeros(1,46); M(1) = 1; tmax = 1000;
% % 
% % K = [
% %     kG1,   %K(1)
% %     kG2,   %K(2)
% %     kG3,   %K(3)
% %     kG4,   %K(4)
% %     kG5,   %K(5)
% %     kp,    %K(6)
% %     kI,    %K(7)
% %     kS,    %K(8)
% %     kC,    %K(9)
% %     kO,    %K(10)
% %     kk1,   %K(11)
% %     kk2,   %K(12)
% %     kk3,   %K(13)
% %     kk4,   %K(14)
% %     kk5 ]; %K(15)
% % 
% % nu = length(K);
% % no_rxns = 10;  % number of reactions (total)
% % h = zeros(no_rxns,1);
% % plot_matrx(1,:) = [t X(1) X(2)]; % plotting matrix
% % tic = 1;  % to construct a matrix for plotting
% % 
% % %% Build function to account for increase in arrestin binding affinity
% % %% with more phosphates bound to melanopsin carboxyl tail.
% % 
% % W = @(n) 1-exp(-n);
% % 
% % %% Build function to account for decrease in G-protein activation
% % %% with more phosphates bound to melanopsin carboxyl tail.
% % 
% % Y = @(n) exp(-n);
% % 
% % while t < tmax
% %     % use cell #'s to construct hazards for the current time step
% %     
% %     % G-Protein Activation
% %     h(1) = Y(0)*K(1)*M(1)*X(1);    % formation of M0*.G.GDP
% %     h(2) = Y(0)*K(2)*M(2);         % dissociation of M0*.G.GDP
% %     h(3) = Y(0)*K(3)*M(2);         % loss of GDP
% %     h(4) = Y(0)*K(4)*M(3);         % gain of GTP
% %     h(5) = Y(0)*K(5)*M(4);         % association of M0*.G with GTP
% %     h(6) = Y(1)*K(1)*M(6)*X(1);
% %     h(7) = Y(1)*K(2)*M(7);
% %     h(8) = Y(1)*K(3)*M(7);
% %     h(9) = Y(1)*K(4)*M(8);
% %     h(10) =Y(1)*K(5)*M(9);
% %     h(11) =Y(2)*K(1)*M(13)*X(1);
% %     h(12) =Y(2)*K(2)*M(14);
% %     h(13) =Y(2)*K(3)*M(14);
% %     h(14) =Y(2)*K(4)*M(15);
% %     h(15) =Y(2)*K(5)*M(16);
% %     h(16) =Y(3)*K(1)*M(20)*X(1);
% %     h(17) =Y(3)*K(2)*M(21);
% %     h(18) =Y(3)*K(3)*M(21);
% %     h(19) =Y(3)*K(4)*M(22);
% %     h(20) =Y(3)*K(5)*M(23);
% %     h(21) =Y(4)*K(1)*M(27)*X(1);
% %     h(22) =Y(4)*K(2)*M(28);
% %     h(23) =Y(4)*K(3)*M(28);
% %     h(24) =Y(4)*K(4)*M(29);
% %     h(25) =Y(4)*K(5)*M(30);
% %     h(26) =Y(5)*K(1)*M(34)*X(1);
% %     h(27) =Y(5)*K(2)*M(35);
% %     h(28) =Y(5)*K(3)*M(35);
% %     h(29) =Y(5)*K(4)*M(36);
% %     h(30) =Y(5)*K(5)*M(37);
% %     h(31) =Y(6)*K(1)*M(41)*X(1);
% %     h(32) =Y(6)*K(2)*M(42);
% %     h(33) =Y(6)*K(3)*M(42);
% %     h(34) =Y(6)*K(4)*M(43);
% %     h(35) =Y(6)*K(5)*M(44);
% %     
% %     % PLC and G-protein activation/inactivation
% %     
% %     h(36) =K(6)*X(3)*X(5);
% %     h(37) =K(7)*X(4)*X(6);
% %     
% %     % Second Messenger Creation
% %     
% %     h(38) =K(8)*X(6)*X(10);
% %     
% %     % Channel Opening
% %     
% %     h(39) =K(9)*X(11)*X(12);
% %     h(40) =K(10)*X(13);
% %     
% %     % Kinase Phosphorylation
% %     
% %     h(41) =K(11)*M(1);
% %     h(42) =
% %     
% %     
% %     
% %     % now turn the hazards into percentages with hw
% %     h_tot = sum(h);
% %     hw = h/h_tot;
% %     
% %     % now pick a random number to step forward in time
% %     %tt = -log(rand(1,1))/h_tot;
% %     tt = exprnd(h_tot);
% %     t = t + tt;
% %     % pick a random number and use the weights to "select" an action
% %     r = rand(1,1);
% %     
% %     % based on r, HSC/P/NK's can increase/decrease/remain
% %     if 0 <= r && r < hw(1) && X(1) > 0
% %         X(1) = X(1) + 1;     %birth of HSC
% %     elseif hw(1) < r && r < hw(1) + hw(2) && X(1) > 0
% %         X(1) = X(1) - 1;     %death of HSC
% %     elseif hw(1) + hw(2) < r && r < hw(1) + hw(2) ...
% %             + hw(3) && X(1) > 0
% %         X(2) = X(2) + 1;     %differentiation into NK
% %     elseif hw(1) + hw(2) + hw(3) < r && r < hw(1) ...
% %             + hw(2) + hw(3) + hw(4) && X(2) > 0 
% %         X(2) = X(2) + 1;     %prolif. of NK
% %     elseif hw(1) + hw(2) + hw(3) + hw(4) < r ...
% %             && r < 1 && X(2) > 0
% %         X(2) = X(2) - 1;     %death of NK
% %     else
% %         X(1) = X(1);   % If the # of cells is zero
% %         X(2) = X(2);   % then nothing can happen
% %     end
% %     
% %     % adjust the plotting matrix with the newfound information
% %     tic = tic + 1;
% %     plot_matrx(tic,:) = [t X(1) X(2)];
% %     
% %     % if the new t < tmax, repeat.  Else plot the results.
% % end
% % 
% % %plotting!
% % hold on
% % plot(plot_matrx(:,1),plot_matrx(:,2),'b:',...
% %     plot_matrx(:,1),plot_matrx(:,3),'r:');
% % axis([0 tmax -1 200]);
% % xlabel('time (/day)'); ylabel('# of cells');



