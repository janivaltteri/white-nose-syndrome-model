## This is a simple example simulation of White Nose Syndrome
## dynamics in bats in a single patch case to accompany
## Lilley, T., Anttila, J., Ruokolainen, L. 2008. Functional Ecology
## Landscape structure and ecology influence the spread of a bat fungal disease

## This can be run directly from command line with:
##  julia wns-simple.jl
## and should generate an output figure (out_inf.png)

## using the package Gaston for plotting
## note: this requires that gnuplot is installed

using Gaston

## define functions

function wakerate(x::Float64, p::Dict{String,Float64})
  y = (x/p["wr_half"])^p["wr_kappa"]
  return y/(1.0 + y)
end

function sleeprate(x::Float64, p::Dict{String,Float64})
  y = (x/p["sr_half"])^p["sr_kappa"]
  return (1.0 - y/(1.0 + y))
end

function infresp(x::Float64, p::Dict{String,Float64})
  y = (x/p["id50"])^p["kappa"]
  return y/(1.0 + y)
end

function outempfromns(nval::Float64)
  mult = 10.0
  add = 10.0
  return mult*nval + add
end

function intempfromns(x::Float64)
  mult = 0.5
  add = 0.0
  return x*mult + add
end

function batgrowth(x::Float64, p::Dict{String,Float64})
  y = (x/p["bg_half"])^p["bg_kappa"]
  return p["rh"]*y/(1.0 + y)
end

function fungalgrowth(x::Float64)
    tmin = 7.0
    tmax = 15.0
    if x < tmin
        return 0.0
    elseif x > tmax
        return 0.0
    else
        maximum = 50.0
        b3 = 0.0377
        c3 = 0.25
        return maximum*(b3*(x-tmin))*(b3*(x-tmin))*(1-exp(c3*(x-tmax)))
    end
end

function hmortality(x::Float64)
    c3 = -0.0001032496
    c2 =  0.0083275607
    c1 = -0.0825981368
    c0 =  0.3468233335
    return c3*x*x*x + c2*x*x + c1*x + c0
end
function stepper(st, par, t)
    nst = zeros(8);
    nval = nvals[t];
    outemp = outempfromns(nval);
    intemp = ct + intempfromns(nval);
    eta1 = sleeprate(outemp,par);
    eta2 = wakerate(outemp,par);
    fgrow = fungalgrowth(intemp) + ct*0.0;
    mgrow = fungalgrowth(intemp) + ct*0.0;
    hmort = hmortality(intemp)*par["myy_h"];
    densdep = (st[1] + st[2] + st[3])/par["kh"];
    dir_tm = par["beta_d"]*st[4]*st[6]/(1.0 + par["eta"]*st[4]);
    nst[1] = batgrowth(outemp,par)*((st[1]+st[2]+st[3])-st[1]*densdep) - par["beta_e"]*st[1]*infresp(st[7],par) + par["delta"]*st[2] - eta1*st[1] + eta2*st[4];
    nst[2] = par["beta_e"]*st[1]*infresp(st[7],par) + par["nyy"]*st[3] - batgrowth(outemp,par)*st[2]*densdep - par["delta"]*st[2] - eta1*st[2] + eta2*st[5];
    nst[3] = eta2*st[6] - eta1*st[3] - par["nyy"]*st[3] - batgrowth(outemp,par)*st[3]*densdep;
    nst[4] = eta1*st[1] - eta2*st[4] - hmort*st[4] - dir_tm;
    nst[5] = eta1*st[2] - eta2*st[5] - hmort*st[5] - par["khi"]*st[5] + dir_tm;
    nst[6] = eta1*st[3] - eta2*st[6] - hmort*st[6] - par["myy_f"]*st[6] + par["khi"]*st[5];
    nst[7] = par["rf"]*st[7]*fgrow - 0.1*par["rf"]*st[7]*(st[7] + par["cfm"]*st[8])/par["kf"] + par["lambda"]*st[2];
    nst[8] = par["rm"]*st[8]*mgrow - 0.1*par["rm"]*st[8]*(st[8] + par["cmf"]*st[8])/par["km"];
    return nst
end

## define variables

pii = 3.14159265;
ct = 9.0;

par = Dict(
    "rh" => 0.0052,
    "kh" => 10.0,
    "rf" => 0.025,
    "kf" => 100.0,
    "beta_e" => 0.05,
    "myy_f" => 0.04,
    "lambda" => 0.5,
    "delta" => 0.075,
    "rm" => 0.025,
    "km" => 100.0,
    "cmf" => 0.0,
    "cfm" => 0.0,
    "myy_h" => 0.01,
    "khi" => 0.01,
    "nyy" => 0.075,
    "beta_d" => 0.001,
    "eta" => 0.2,
    "id50" => 140.0,
    "kappa" => 5.0,
    "wr_half" => 13.0,
    "wr_kappa" => 20.0,
    "sr_half" => 7.0,
    "sr_kappa" => 10.0,
    "bg_half" => 10.0,
    "bg_kappa" => 6.0);

nstates = 8;
st = zeros(nstates);
ch = zeros(nstates);
st[1] = 2.0;
st[2] = 0.0;
st[3] = 0.0;
st[7] = 100.0;

steps = 10*3650;
res = 10;
dt = 0.1;
times = collect(0.0:dt:(steps*dt - dt))

nvals = [cos(2.0*pii*t/365.0)*0.4 for t in times]

out = zeros(floor(steps/res),nstates);
out[1,1] = st[1];
out[1,2] = st[2];
out[1,3] = st[3];
out[1,7] = st[7];

meta = zeros(floor(steps/res),2);
meta[1,1] = times[1];
meta[1,2] = nvals[1];

## run simulation loop

counter = 2
for t = 1:(steps-1)
    ch = stepper(st,par,t)
    for i = 1:nstates
        st[i] = st[i] + dt*(ch[i]);
    end
    if mod(t,res) == 0
        for i = 1:nstates
            out[counter,i] = st[i];
        end
        meta[counter,1] = times[t];
        meta[counter,2] = nvals[t];
        counter += 1;
    end
end

plot(meta[:,1],out[:,1])
plot!(meta[:,1],out[:,2],color="red")
plot!(meta[:,1],out[:,3],color="green")
plot!(meta[:,1],out[:,4])
plot!(meta[:,1],out[:,5],color="red")
plot!(meta[:,1],out[:,6],color="green")

printfigure(term="png",outputfile="out_inf.png")
