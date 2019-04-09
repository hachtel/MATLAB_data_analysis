function[coal, coal_err] = self_HOM_coal(tau, numer_g2_para, numer_g2_perp, window, back, second_pulse)
%This function calculates coalescence from a self-HOM measurement.
%We assume the taus are the same for both measurements and that g(2)(0) = 0

%window is how wide you want to sum over
%back is location used to perform background subtraction
%second_pulse is the location of the pulse used for normalization

%Find indices for a window around 0
tau_low = find(tau >=-window, 1);
tau_high = find(tau >=window, 1);

%Find relevant indices for background
tau_low_back = find(tau >=(back-window), 1);
tau_high_back = find(tau >=(back+window), 1);

%find indices for second pulse
tau_2_low = find(tau >= (second_pulse - window), 1);
tau_2_high = find(tau >=(second_pulse + window), 1);

%Now perform the summations
para = sum(numer_g2_para(1,tau_low:tau_high));
perp = sum(numer_g2_perp(1,tau_low:tau_high));
back_para = sum(numer_g2_para(1,tau_low_back:tau_high_back));
back_perp = sum(numer_g2_perp(1,tau_low_back:tau_high_back));
norm_para = sum(numer_g2_para(1,tau_2_low:tau_2_high));
norm_perp = sum(numer_g2_perp(1,tau_2_low:tau_2_high));

%To calculate the error it is easier to redefine these parameters
a = para;
b = back_para;
c = norm_para;
d = perp;
e = back_perp;
f = norm_perp;

%error: we did coal = ((a-b)/(c-b))/(d-e)/(f-e)) = (a-b)/(c-b)*(f-e)/(d-e)
%For ease of reading, define the partial derivates as di where i is the
%variable of interest, corresponding to the list immediately above.
da = (-e + f)/((-b + c)*(d - e));
db = ((a - b)*(-e + f))/((-b + c)^2*(d - e)) - (-e + f)/((-b + c)*  (d - e));
dc = -(((a - b)*(-e + f))/((-b + c)^2*(d - e)));
dd = -(((a - b)*(-e + f))/((-b + c)*(d - e)^2));
de = -((a - b)/((-b + c)*(d - e))) + ((a - b)*(-e + f))/((-b + c)*(d - e)^2);
df = (a - b)/((-b + c)*(d - e));
%Calculate the coalescence
coal = 1-((a - b)/(c - b))/((d - e)/(f - e));
%And the coalescence error
coal_err = (1-((a - b)/(c - b))/((d - e)/(f - e)))*sqrt(da^2*a + db^2*b + dc^2*c + dd^2*d + de^2*e + df^2*f);


