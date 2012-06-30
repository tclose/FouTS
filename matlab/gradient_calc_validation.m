function gradient_calc_validation

%----------------------------------------------------%
%----------------------------------------------------%
%                   Constants                        %
%----------------------------------------------------%
%----------------------------------------------------%

%Number of 'time' (free parameter of the strand parameterisation) samples.
T = 20

%Degree of Fourier Descriptors
D = 3

%DW sample orientations
U = [ 0.9413   -0.0407   -0.3351;
   -0.8846   -0.3834   -0.2654;
   -0.1714   -0.7779    0.6046;
    0.4040   -0.2802   -0.8708;
    0.1593   -0.7891    0.5932;
    0.1960   -0.3076    0.9311;
   -0.0815    0.2989    0.9508;
    0.2384    0.7604    0.6042;
    0.6346    0.7660   -0.1025;
    0.9924   -0.1223   -0.0151;
    0.2629   -0.5636   -0.7831;
    0.1511    0.9332    0.3260;
   -0.4723   -0.7962   -0.3783;
    0.3909         0    0.9204;
    0.6807    0.5251    0.5108;
    0.2602   -0.9654   -0.0149;
    0.5118   -0.2842    0.8107;
    0.6354    0.4318   -0.6402;
    0.5055   -0.8350   -0.2174;
    0.3404   -0.5611    0.7545;
   -0.3686   -0.8923    0.2605;
   -0.4293   -0.5237   -0.7358;
    0.7402   -0.0887    0.6665;
   -0.7366   -0.5749    0.3561;
   -0.6702    0.4612   -0.5815;
   -0.8144   -0.2115   -0.5404;
    0.0739    0.9973         0;
   -0.1919    0.9231    0.3332;
   -0.7180   -0.6653   -0.2047;
    0.2784   -0.9094    0.3090;
   -0.5392    0.8337   -0.1189;
   -0.5836   -0.2303   -0.7787;
   -0.0742   -0.5672   -0.8203;
   -0.4803    0.7237   -0.4955;
   -0.0686         0   -0.9976;
    0.6979   -0.6031   -0.3863;
   -0.0887    0.7850    0.6131;
    0.6301   -0.0634   -0.7739;
   -0.2544   -0.3018   -0.9188;
   -0.9045   -0.2937    0.3093;
    0.9044   -0.4125    0.1089;
    0.4737    0.2298   -0.8502;
    0.2616    0.0001   -0.9652;
    0.3362    0.5345   -0.7754;
    0.9786    0.2057         0;
    0.1428    0.3091   -0.9402;
    0.5723   -0.4776   -0.6666;
   -0.0492   -0.9405    0.3361;
    0.9566    0.0391    0.2889;
    0.7269   -0.6210    0.2932;
   -0.0076    0.5804   -0.8143;
    0.8580    0.5109   -0.0530;
   -0.3915   -0.9169   -0.0770;
    0.4133   -0.7417   -0.5283;
    0.7995    0.1435   -0.5833;
    0.8747   -0.2627    0.4072;
   -0.4816   -0.7070    0.5179;
   -0.7969    0.2592    0.5457;
    0.7540   -0.6540   -0.0621;
   -0.8983    0.3722    0.2335];

  %Pregenerated Fourier Descriptor basis set;
  PSI = fourier_descr_matrix(T, D, 0);

  %Pregenerated derived Fourier Descriptor w.r.t 'time' basis set
  dPSI = d_fourier_descr_matrix(T, D, 0);
     
  %Signal response model coefficients
%   Resp = [0.3475; -1.6274; 3.2524; -3.0952; 1.1302];  
  
  Resp = [0.264387; -1.50073; 1.89986; -1.28088; 0.366885];
    
  %Voxel centre
  VOX = [0.0, 0.0, 0.0];  
  
  %Gaussian blur standard deviation
  ALPHA = 1.0;
  
  
  
%----------------------------------------------------%
%----------------------------------------------------%
%                Default Indices                     %
%----------------------------------------------------%
%----------------------------------------------------%
  

%   %Response coefficient index
   j = 4;
 
%   %Strand free parameter, 'time', index.
   t = 10
% 
%   %Strand degree index  
   m = 2  
% 
%   %Strand spatial dimension index
   d = 1
% 
%   %DW Sample orientation index
   i = 1
% 
%   %Voxel index
   n = 1
  
  
%----------------------------------------------------%
%----------------------------------------------------%
%               Primary Functions                    %  
%----------------------------------------------------%
%----------------------------------------------------%  

  %Point 't', along strand path.
  function lambda_t = Lambda(V,t)

    lambda = PSI * V;
    
    lambda_t = lambda(t,:)';
    
  end

  %Gradient w.r.t V(m,d)
  function dlambda_t = dLambda(V,t,m)

    dlambda_t = PSI(t,m);
    
  end


  %Tangent at point 't' along strand path 
  function glambda_t = GLambda(V,t)

    glambda = dPSI * V;
    
    glambda_t = glambda(t,:)';
    
  end

  %Gradient w.r.t V(m,d)
  function dglambda_t = dGLambda(V,t,m)

    dglambda_t = dPSI(t,m);
    
  end



  %Signal response model coefficient 'j'
  function r = R(j)
    
    r = Resp(j+1);
    
  end

  %Voxel n
  function x = X(n)
    
    x = VOX(n,:)';
    
  end


  %Gaussian smoothing kernel normalising scalar.
  function scalar = ConvScalar
    
    scalar = 1/( sqrt(2) * pi * ALPHA );
    
  end

  
%----------------------------------------------------%
%----------------------------------------------------%
%              Secondary Functions                   %  
%----------------------------------------------------%
%----------------------------------------------------%  


  % %tang_dw_dot% Dot product between strand tangent and DW-sample orientation 'i'.
  function rho_ti = Rho(V, t, i)
    
   rho_ti = U(i,:) * GLambda(V,t);
      
  end


  %Gradient w.r.t V(m,d)
  function drho_ti = dRho(V, t, i, d, m)
    
   drho_ti = U(i,d) * dGLambda(V,t,m);
      
  end


  % %param_speed2% (Parameterisation speed)^2 at point 't' along strand path
  function zeta_t = Zeta(V, t)
    
    zeta_t = GLambda(V,t)' * GLambda(V,t);

  end

  %Gradient w.r.t V(m,d)
  function dzeta_t = dZeta(V, t, d, m)
    
    gLambda_Vt = GLambda(V,t);
    
    dzeta_t = 2 * dGLambda(V,t,m) * gLambda_Vt(d);

  end

  % %strand_to_vox_disp% Displacement along dimension 'd', between voxel 'n' and point 't' along strand path, 
  %normalised by the standard deviation of the Gaussian smoothing kernel. 
  function delta_t = Delta(V, t, n, d)
    
    lambda_Vt = Lambda(V, t);
    Xn = X(n);
    delta_t = (lambda_Vt(d) - Xn(d))/ALPHA;
    
  end

  %Gradient w.r.t V(m,d)
  function ddelta_t = dDelta(V, t, n, d, m)
    
    ddelta_t = dLambda(V,t,m)/ALPHA;
    
  end


  
%----------------------------------------------------%
%----------------------------------------------------%
%               Tertiary Functions                   %  
%----------------------------------------------------%
%----------------------------------------------------%




  %Dot product between strand tangent and DW-sample orientation 'i'.
  function qrho_ti = QRho(V, t, i)
    
   qrho_ti = Rho(V,t,i) ^ (2*j);
      
  end


  %Gradient w.r.t V(m,d)
  function dqrho_ti = dQRho(V, t, i, d, m)
       
   dqrho_ti = dRho(V,t,i,d,m) * (2*j) * Rho(V,t,i) ^ (2*j-1);
      
  end


  %(Parameterisation speed)^2 at point 't' along strand path
  function qzeta_t = QZeta(V, t)
    
    qzeta_t = Zeta(V,t) ^ ( (1-2*j)/2 );

  end

  %Gradient w.r.t V(m,d)
  function dqzeta_t = dQZeta(V, t, d, m)
    
    dqzeta_t = dZeta(V,t,d,m) * ( (1-2*j)/2 ) * Zeta(V,t) ^ ( (-1-2*j)/2 );

  end

  %Displacement along dimension 'd', between voxel 'n' and point 't' along strand path, 
  %normalised by the standard deviation of the Gaussian smoothing kernel. 
  function qdelta_t = QDelta(V, t, n, d)
    
    qdelta_t = exp(-Delta(V,t,n,d)^2);
    
  end

  %Gradient w.r.t V(m,d)
  function dqdelta_t = dQDelta(V, t, n, d, m)
    
    dqdelta_t = -2 * dDelta(V,t,n,d,m) * Delta(V,t,n,d) * exp (-Delta(V,t,n,d)^2);
    
  end




%----------------------------------------------------%
%----------------------------------------------------%
%                Partial Functions                   %  
%----------------------------------------------------%
%----------------------------------------------------%  


  %Signal component due to the 'j'th response harmonic from a point 't' along the strand that shows up in 
  % the 'i'th DW-sample orientation, at a voxel 'n', displaced along dimension 'd' 
  function q = Q(V,t,j,i,n,d)
    
    q = R(j) * ConvScalar * QRho(V,t,i) * QZeta(V,t) * QDelta(V,t,n,d);
    
  end

  %Gradient w.r.t V(m,d)
  function q = dQ(V,t,j,i,n,d,m)
    
    q = R(j) * ConvScalar * (dQRho(V,t,i,d,m) * QZeta(V,t) * QDelta(V,t,n,d) + QRho(V,t,i) * dQZeta(V,t,d,m) * QDelta(V,t,n,d) + QRho(V,t,i) * QZeta(V,t) * dQDelta(V,t,n,d,m));
    
  end


  %Signal component due to the 'j'th response harmonic from a point 't' along the strand that shows up in 
  % the 'i'th DW-sample orientation, at a voxel 'n', displaced along dimension 'd' 
  function q = seg_signal(V,t,j,i)
    
    
    q = R(j) * QRho(V,t,i) * QZeta(V,t);
    
  end

  %Gradient w.r.t V(m,d)
  function q = seq_gradient(V,t,j,i,n,d,m)
    
    q = dQRho(V,t,i,d,m) ;%* QZeta(V,t);%R(j) *  ... + QRho(V,t,i) * dQZeta(V,t,d,m));
    
  end



%----------------------------------------------------%
%----------------------------------------------------%
%                Partial Functions                   %  
%----------------------------------------------------%
%----------------------------------------------------%  




%   function q_complete = Q_complete(V,i,n)
%     
%     q_complete = 0;
%     
%     for t = 1:T
%       
%       for j = 0:1:(size(Resp,1)-1)
%         
%         for d = 1:3
%         
%           q_complete = q_complete + Q(V,t,j,i,n,d);
%         
%         end  
%           
%       end
%       
%     end    
%     
%     
%   end
% 
% 
%   function dq_complete = dQ_complete(V,i,n)
%     
%     dq_complete = zeros(size(V));
%     
%     for m = 1:size(V,1)
%       
%       for d = 1:3      
% 
%         for t = 1:T
% 
%           for j = 0:1:(size(Resp,1)-1)
% 
%               dq_complete(m,d) = dq_complete(m,d) + dQ(V,t,j,i,n,d,m);
% 
%           end  
% 
%         end
% 
%       end   
%     
%     end
%     
%   end

%----------------------------------------------------%
%----------------------------------------------------%
%                   Validation                       %  
%----------------------------------------------------%
%----------------------------------------------------%  

%   %small increment
   h = 0.000001;

  
  %Strand Fourier Descriptor (Chung, 2008) coefficients
  V_INIT = [  -0.00957484, -0.00342075, -0.02118;...
              -0.209237, 0.0614428, -0.0497531;...
               0.00225339, 0.00672736, 0.0199366 ]
    

    
  M = size(V_INIT,1)
%   
%   Q_complete_ = Q_complete(V_INIT,i,n)
%   
%   dQ_complete_ref = zeros(M,3);
%   
%   for m = 1:M
%     for d = 1:3
%    
%       Vh = V_INIT;
%       Vh(m,d) = Vh(m,d) + h;
%       
%       dQ_complete_ref(m,d) = (Q_complete(Vh,i,n) - Q_complete(V_INIT,i,n))/h;
%       
%     end
%   end
%   
%   dQ_complete_ref
%   
%   dQ_complete_  = dQ_complete(V_INIT,i,n)
 
  

  H = zeros(D,3);
  H(m,d) = h;
  Vh = V_INIT + H

  
disp('-------------- Response coeff ------------');       
  
  R(j)
  
disp('-------------- Segment -----------------');       

  Lambda_Vt = Lambda(V_INIT,t)

  GLambda_Vt = GLambda(V_INIT,t)


  
disp('------------- Complete Function  ------------------');       
  
%  Evaluation of the signal component.
  Q_Vtjind      = Q(V_INIT,t,j,i,n,d)

%   Evalutation of the estimated gradient.
  dQ_Vtjindm_ref = ( Q(Vh,t,j,i,n,d) - Q(V_INIT,t,j,i,n,d))/h
   
  
%   Evaluation of the calculated gradient.
   dQ_Vtjindm = dQ(V_INIT,t,j,i,n,d,m)  

disp('--------------- Rho  ----------------');     

%   %Evaluation of the calculated gradient.
   Rho_Vti = Rho(V_INIT,t,i)


%   %Evalutation of the estimated gradient.
%    dRho_Vtindm_ref = ( Rho(Vh,t,i) - Rho(V_INIT,t,i))/h
 
%   %Evaluation of the calculated gradient.
   dRho_Vtindm = dRho(V_INIT,t,i,d,m)


disp('--------------- Zeta ----------------');  

  Zeta_Vt = Zeta(V_INIT,t)

%   %Evalutation of the estimated gradient.
%    dZeta_Vt_ref = ( Zeta(Vh,t) - Zeta(V_INIT,t))/h
% 
%   %Evaluation of the calculated gradient.
   dZeta_Vtdm = dZeta(V_INIT,t,d,m)
%

disp('--------------- Delta ----------------');  

%   %Evalutation of the estimated gradient.
%    dDelta_Vtnd_ref = ( Delta(Vh,t,n,d) - Delta(V_INIT,t,n,d))/h
% 
%   %Evaluation of the calculated gradient.
   dDelta_Vtndm = dDelta(V_INIT,t,n,d,m)  
   
disp('--------------- Q Rho ----------------');

   QRho_Vti = QRho(V_INIT, t,i)

%   %Evaluation of the calculated gradient.
   dQRho_Vtindm = dQRho(V_INIT,t,i,d,m)
   
disp('-------------- QZeta -----------------');     
   
  QZeta_Vt = QZeta(V_INIT,t)

%   %Evaluation of the calculated gradient.
   dQZeta_Vtdm = dQZeta(V_INIT,t,d,m)   
   

disp('-------------- Segment Gradient (minus interpolation)-----------');         
  
  seg_signal_Vtji = seg_signal(V_INIT,t,j,i)
    
  seg_gradient_Vtjindm = seq_gradient(V_INIT,t,j,i,n,d,m)  
  
  seg_gradient_Vtjindm_ref = (seg_signal(Vh,t,j,i) - seg_signal(V_INIT,t,j,i))/h  
     
   
end

