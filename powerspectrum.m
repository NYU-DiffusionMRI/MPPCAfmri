
clear

p = phantom(60);
pr = repmat(p,[1,1,3,100]);
pn = pr + randn(size(pr));
pn_ = reshape(pn, [60*60*3,100]);
[pdn_, s_] = MP(pn_);
sigma = s_.*ones(60,60,3);
pdn = reshape(pdn_,[60,60,3,100]);

imorig = pn;
imdn = pdn;

% strongly reccomendd a mask
[ps, w, ac, t] = powerspectrum_greg(imorig, imdn, sigma);

figure;
subplot(1,2,1)
plot(w, mean(ps,2))
subplot(1,2,2)
plot(t, mean(ac,2));

function [powerspectrum,w,autocorrelation,t] = powerspectrum_greg(original, denoised, sigma)
    [sx, sy, sz, T] = size(original);
    original_ = reshape(original, [T, sx*sy*sz]);
    denoised_ = reshape(denoised, [T, sx*sy*sz]);
    sigma_ = reshape(sigma, [1, sx*sy*sz]);
    
    dt = 1; 
    t = -T/2:dt:T/2-dt; 
    dw = 2*pi/T; 
    w = -pi/dt:dw:pi/dt-dw;
    
    ep = (denoised_ - original_)./(sigma_ + eps);
    
    epf = fftshift(fft(ep))*dt; 
    Gamma = conj(epf).*epf / T; 
    Gamma_t = ifftshift(ifft(Gamma))/dt;
    
    w = w/2/pi;
    powerspectrum = Gamma;
    autocorrelation = Gamma_t;
end    

