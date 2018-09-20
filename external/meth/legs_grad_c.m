function gradbasis=legs_grad_c(x,dir,n,scale)
% usage: gradbasis=legs_grad_c(x,dir,n,scale)
%        
%         only Y_nm/r^(2n+1)
% returns the values and directional derivatives  of (n+1)^2-1 basis functions 
% constructed from spherical harmonics at locations given in x and, for the 
% gradients, for (in general non-normalized) directions given in dir.   
% 
% input: x      set of N locations given as an Nx3 matrix 
%        dir    set of N direction vectors given as an Nx3 matrix 
%                  (dir is not normalized (it hence can be a dipole moment))
%        n       order of spherical harmonics 
%
% output: basis: Nx((n+1)^2-1)  matrix containing in the j.th  row the real 
%                and imaginary parts of r^kY_{kl}(theta,Phi)/(N_{kl}*scale^k) ( (r,theta,phi) 
%                are the spherical coordinates corresponding to  the j.th row in x) 
%                for k=1 to n and l=0 to k 
%                the order is:
%                          real parts for k=1 and l=0,1 (2 terms) 
%                  then    imaginary parts for k=1 and l=1 (1 term) 
%                  then    real parts for k=2 and l=0,1,2 (3 terms) 
%                  then    imaginary parts for k=2 and l=1,2 (2 term) 
%                              etc.
%                   the spherical harmonics are normalized with
%                   N_{kl}=sqrt(4pi (k+l)!/((k-l)!(2k+1)))
%                    the phase does not contain the usual (-1)^l term !!! 
%                   scale is constant preferably set to the avererage radius                   
%
%         gradbasis: Nx((n+1)^2-1) matrix containing in the j.th row the scalar 
%                     product of the gradient of the former with the j.th row of dir
%             
% CC Guido Nolte 
%

[n1,n2]=size(x);

comi=sqrt(-1);

normalize=ones(n,n+1);
facto=ones(2*n+2,1);for i=3:2*n+2,facto(i)=facto(i-1)*(i-1);end;
for i=1:n
    for j=1:i+1
        if j==1
            normalize(i,j)=scale^i*i*sqrt(2*facto(i+j-1+1)/facto(i-j+1+1)/(2*i+1));
        else 
            normalize(i,j)=scale^i*i*sqrt(facto(i+j-1+1)/facto(i-j+1+1)/(2*i+1));
        end
    end
end
normalize=reshape(repmat(reshape(normalize,n*(n+1),1),1,n1),n,n+1,n1);



normalize_b=ones(n,n+1);
for i=1:n
    for j=1:i+1
            normalize_b(i,j)=scale^(2*i+1);
    end
end
normalize_b=reshape(repmat(reshape(normalize_b,n*(n+1),1),1,n1),n,n+1,n1);



rad=sqrt(x(:,1).^2+x(:,2).^2+x(:,3).^2);
phi=angle(x(:,1)+comi*x(:,2));
costheta=x(:,3)./(rad+eps);

ms=0:n;
ns=1:n;

shiftfactors=zeros(n,n+1);
shiftminusfactors=zeros(n,n+1);
for in=1:n; for im=1:n+1
        shiftfactors(in,im)=in+im-1;
        shiftminusfactors(in,im)=(in+im-2)*(in+im-1);
    end;end;
for in=1:n
    shiftminusfactors(in,1)=1;
end;


leg0=zeros(n,n+1,n1);

for i=1:n
 p=legendre(i,costheta);
 leg0(i,1:i+1,:)=p;
end
 leg0=leg0.*reshape(repmat((-1).^ms,n,n1),n,n+1,n1);
ephi=exp(comi*ms'*phi');

leg0=leg0.*reshape(repmat(reshape(ephi,1,(n+1)*n1),n,1),n,n+1,n1)...
         .*reshape(repmat((repmat(rad,1,n).^repmat(ns,n1,1))',n+1,1),n,n+1,n1);

     onesx=zeros(1,n+1,n1);
     onesx(1,1,:)=ones(1,1,n1);
     
legshift=[onesx;leg0(1:n-1,:,:)];
legshiftminus=[-conj(legshift(:,2,:)),legshift(:,1:n,:)].*reshape(repmat(reshape(shiftminusfactors,n*(n+1),1),n1,1),n,n+1,n1);
legshiftplus=-[legshift(:,2:n+1,:),zeros(n,1,n1)];
legshift=legshift.*reshape(repmat(reshape(shiftfactors,n*(n+1),1),n1,1),n,n+1,n1);

dirp=[(dir(:,1)+dir(:,2)/comi)/2,(dir(:,1)-dir(:,2)/comi)/2,dir(:,3)];



gradleg= reshape(repmat(transpose(dirp(:,1)),n*(n+1),1),n,n+1,n1).*legshiftplus...
        +reshape(repmat(transpose(dirp(:,2)),n*(n+1),1),n,n+1,n1).*legshiftminus...
        +reshape(repmat(transpose(dirp(:,3)),n*(n+1),1),n,n+1,n1).*legshift;
    
leg0=leg0./normalize;
gradleg=gradleg./normalize;

 xdotdir_ori= (x(:,1).*dir(:,1)+x(:,2).*dir(:,2)+x(:,3).*dir(:,3))';
 xdotdir= reshape(repmat(xdotdir_ori,n*(n+1),1),n,n+1,n1);
      

leg0_b=leg0./reshape(repmat((repmat(rad,1,n).^repmat(2*ns+1,n1,1))',n+1,1),n,n+1,n1);
gradleg_b= gradleg./reshape(repmat((repmat(rad,1,n).^repmat(2*ns+1,n1,1))',n+1,1),n,n+1,n1)...
             -leg0./reshape(repmat((repmat(rad,1,n).^repmat(2*ns+3,n1,1))',n+1,1),n,n+1,n1)...
                  .*reshape(repmat((repmat(2*ns+1,n1,1))',n+1,1),n,n+1,n1).*xdotdir;
% 
% rad=rad
% scale=scale
%          
% leg0_b=leg0_b.*normalize_b;
gradleg_b=gradleg_b.*normalize_b;

 % 
%          x_a=reshape(gradleg,n*(n+1),n1);
%          x_b=reshape(gradleg_b,n*(n+1),n1);
%          ta=norm(x_a,'fro')
%          tb=norm(x_b,'fro')
%      
           
    
%     basis=zeros(n1,2*(n+1)^2-1);
 gradbasis=zeros(n1,(n+1)^2);

 %  for i=1:n
% %      basis(:,i^2:i^2+i)=(reshape(real(leg0(i,1:i+1,:)),i+1,n1))';
% %      basis(:,i^2+i+1:(i+1)^2-1)=(reshape(imag(leg0(i,2:i+1,:)),i,n1))';
%      gradbasis(:,i^2:i^2+i)=(reshape(real(gradleg(i,1:i+1,:)),i+1,n1))';
%      gradbasis(:,i^2+i+1:(i+1)^2-1)=(reshape(imag(gradleg(i,2:i+1,:)),i,n1))';    
%  end
    
%     basis(:,n_old+1)=1./rad*scale;
    gradbasis(:,1)=-1./rad.^3.*xdotdir_ori'*scale;
    
    k=1;
    
 for i=1:n
%      basis(:,k+i^2:k+i^2+i)=(reshape(real(leg0_b(i,1:i+1,:)),i+1,n1))';
%      basis(:,k+i^2+i+1:k+(i+1)^2-1)=(reshape(imag(leg0_b(i,2:i+1,:)),i,n1))';
     gradbasis(:,k+i^2:k+i^2+i)=(reshape(real(gradleg_b(i,1:i+1,:)),i+1,n1))';
     gradbasis(:,k+i^2+i+1:k+(i+1)^2-1)=(reshape(imag(gradleg_b(i,2:i+1,:)),i,n1))';    
 end
    


return;
