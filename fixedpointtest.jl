### A Pluto.jl notebook ###
# v0.14.7

using Markdown
using InteractiveUtils

# ╔═╡ 8fdf5004-b7c5-4575-b455-2ec40748f355
using Random

# ╔═╡ 927a8e03-2db7-4d35-b391-5d5cd13f6ea4
using LinearAlgebra

# ╔═╡ 1764fb52-c62d-11eb-0b88-9f632da0c4f7
"""
return a list of lists from the infilename cof data file
    dyear,
    n,m,G,H,SecG,SecH,
    ...

    Args:
        infilename(string ending in .COF): the .COF file that contains the
            WMM coefficents, download this from https://www.ngdc.noaa.gov/geomag/WMM/DoDWMM.shtml
        maxdegree(positive integer): maximum degree
"""
function parseescof(infilename, maxdegree)
	s=""
	open(infilename) do f
		s= read(f,String)
	end
	dyear= parse(Float64,split(s)[1])
	n= map(x->parse(Int,x),split(s)[4:6:end-2])
	m= map(x->parse(Int,x),split(s)[5:6:end-2])
	G= map(x->parse(Float64,x),split(s)[6:6:end-2])
	H= map(x->parse(Float64,x),split(s)[7:6:end-2])
	SecG= map(x->parse(Float64,x),split(s)[8:6:end-2])
	SecH= map(x->parse(Float64,x),split(s)[9:6:end-2])
	return (dyear= dyear, n=n, m=m, G=G, H=H, SecG=SecG, SecH=SecH) 
end

# ╔═╡ 57f856ab-c04e-46a6-883d-44eb39c248fe


# ╔═╡ d2647b9f-5c7d-4f7c-a559-204ad8081ad3
N= 12

# ╔═╡ 22dedee1-7710-41de-82d8-453ffd9c798a
L= ((N+1) * N)÷2 + N

# ╔═╡ 2ec9b6be-2b88-4360-9585-893e01a8aa52
2^15

# ╔═╡ 0a0477b3-ff5f-4c13-9168-57810ee2549f
R= 6371200.0

# ╔═╡ 85e7c3fc-0dd8-4b9d-a224-f6a97db2f210
const cofs=parseescof("test_codegen/WMM2015.cof",12)

# ╔═╡ 5eb05e2f-08f6-49c7-aede-4bc8a4f885e4
string(round.(Int,cofs.H*10))

# ╔═╡ db18d18c-1937-457c-a18a-c55a5f6667db
string(round.(Int,cofs.SecH*10))

# ╔═╡ f60ba0ac-aab6-4fe8-a281-d50440cf09bf
function C(n,m)
	i= ((n+1) * n)÷2 + m
	return cofs.G[i]
end

# ╔═╡ 8c5baa99-7043-4658-b1bf-c4eddde7a734
function S(n,m)
	i= ((n+1) * n)÷2 + m
	return cofs.H[i]
end

# ╔═╡ 3c8e550d-e407-499c-a48a-d4a31351c1be
S(12,12)

# ╔═╡ b0a5031b-6890-42eb-9655-52d57eddc1c0
((N+1) * N)÷2 + N

# ╔═╡ 16f5d989-edf6-4c18-aa95-963ed934a33a
length(cofs.G)

# ╔═╡ 63cf4e77-cda9-442f-b78d-f8feadadfe81
begin
    x= 1111164.8708100126;
    y= 0.0;
    z= 6259542.961028692;
    truthx= -15978.489161206862;
    truthy= -445.9;
    truthz= -52454.5672160699;
end

# ╔═╡ bfa7f0f0-60c1-49de-b337-252548ee34e3
-1.5978489161206863e-05*1E9

# ╔═╡ e24a3ce2-93d6-4e2b-85fb-008edc98ab31
function fullmagfield(x,y,z,N,R,C= C,S= S)
    #step 1 get Vs and Ws
    px=0.0
    py=0.0
    pz=0.0
    r= sqrt(x*x+y*y+z*z)
    r2= x*x+y*y+z*z
    Vtop= R/r
    Wtop= 0.0
    Vprev= 0.0
    Wprev= 0.0
    Vnm= Vtop
    Wnm= Wtop
    a= z*R/r2
    b= R*R/r2
    c= x*R/r2
    d= y*R/r2
    
    for m in 0:N+1
        sqrtmult= (1/sqrt(1))*(1/sqrt(2*m+1))
        lastsqrtmult= 0
        for n in m:N+1
            if(n==m)
                if(m!=0)
                    if(m==1)
                        nm_mult=1
                    else
                        nm_mult= (2*m-1)*(1/sqrt(2*m-1))*(1/sqrt(2*m))
					end
                    temp= Vtop
                    Vtop= nm_mult*(c*Vtop - d*Wtop)
                    Wtop= nm_mult*(c*Wtop + d*temp) 
                    Vprev= 0.0
                    Wprev= 0.0
                    Vnm= Vtop
                    Wnm= Wtop
				end
            else
                temp= Vnm
                Vnm= (2*n-1)*sqrtmult*a*Vnm - (n+m-1)*(n-m-1)*sqrtmult*lastsqrtmult*b*Vprev
                Vprev= temp
                temp= Wnm
                Wnm= (2*n-1)*sqrtmult*a*Wnm - (n+m-1)*(n-m-1)*sqrtmult*lastsqrtmult*b*Wprev
                Wprev= temp
                lastsqrtmult= sqrtmult
                sqrtmult= (1/sqrt(n+1-m))*(1/sqrt(n+1+m))
			end
            if (m<N && n>=m+2)
                f= 1
                if m==0
                    f= sqrt(2)
				end
                px+= 0.5*f*(C(n-1,m+1)*Vnm+S(n-1,m+1)*Wnm)*sqrt((n-m)*(n-m-1))
                py+= 0.5*f*(-C(n-1,m+1)*Wnm+S(n-1,m+1)*Vnm)*sqrt((n-m)*(n-m-1))
			end
            if (n>=2 && m>=2)
                px+= 0.5*(-C(n-1,m-1)*Vnm-S(n-1,m-1)*Wnm)*sqrt((n+m)*(n+m-1))
                py+= 0.5*(-C(n-1,m-1)*Wnm+S(n-1,m-1)*Vnm)*sqrt((n+m)*(n+m-1))
			end
            if (m==1 && n>=2)
                px+= -C(n-1,0)*Vnm*sqrt((n+1)*(n)/2)
                py+= -C(n-1,0)*Wnm*sqrt((n+1)*(n)/2)
			end
            if (n>=2 && n>m)
                pz+= (-C(n-1,m)*Vnm-S(n-1,m)*Wnm)*sqrt((n-m)*(m+n))
			end
		end
	end    
    return (px,py,pz)
end

# ╔═╡ 048c58bc-b114-40e1-9506-cac194941cc1
fullmagfield(1111164.8708100126,0.0,6259542.961028692,N,R,C,S)

# ╔═╡ 1320035e-4ffd-45c1-a502-871fc03bbf19
fullmagfield(1128529.6885767058,0.0,6358023.736329913,N,R,C,S)

# ╔═╡ 4c60844a-bba9-4192-870f-a2a6c80c6f99
"""
Coeffient to feild vectors, This times the coeffients [G,H] give the feild.
"""
function fieldcoeffs(x,y,z,N,R)
    #step 1 get Vs and Ws
	i(n,m)= ((n+1) * n)÷2 + m
	L= i(N,N)
	cofsx= zeros(2L)
	cofsy= zeros(2L)
	cofsz= zeros(2L)
    r= sqrt(x*x+y*y+z*z)
    r2= x*x+y*y+z*z
    Vtop= R/r
    Wtop= 0.0
    Vprev= 0.0
    Wprev= 0.0
    Vnm= Vtop
    Wnm= Wtop
    a= z*R/r2
    b= R*R/r2
    c= x*R/r2
    d= y*R/r2
    
    for m in 0:N+1
        sqrtmult= (1/sqrt(1))*(1/sqrt(2*m+1))
        lastsqrtmult= 0
        for n in m:N+1
            if(n==m)
                if(m!=0)
                    if(m==1)
                        nm_mult=1
                    else
                        nm_mult= (2*m-1)*(1/sqrt(2*m-1))*(1/sqrt(2*m))
					end
                    temp= Vtop
                    Vtop= nm_mult*(c*Vtop - d*Wtop)
                    Wtop= nm_mult*(c*Wtop + d*temp) 
                    Vprev= 0.0
                    Wprev= 0.0
                    Vnm= Vtop
                    Wnm= Wtop
				end
            else
                temp= Vnm
                Vnm= (2*n-1)*sqrtmult*a*Vnm - (n+m-1)*(n-m-1)*sqrtmult*lastsqrtmult*b*Vprev
                Vprev= temp
                temp= Wnm
                Wnm= (2*n-1)*sqrtmult*a*Wnm - (n+m-1)*(n-m-1)*sqrtmult*lastsqrtmult*b*Wprev
                Wprev= temp
                lastsqrtmult= sqrtmult
                sqrtmult= (1/sqrt(n+1-m))*(1/sqrt(n+1+m))
			end
            if (m<N && n>=m+2)
                f= 1
                if m==0
                    f= sqrt(2)
				end
				cofsx[i(n-1,m+1)]+= 0.5*f*Vnm*sqrt((n-m)*(n-m-1))
				cofsx[i(n-1,m+1)+L]+= 0.5*f*Wnm*sqrt((n-m)*(n-m-1))
				cofsy[i(n-1,m+1)]+= -0.5*f*Wnm*sqrt((n-m)*(n-m-1))
				cofsy[i(n-1,m+1)+L]+= 0.5*f*Vnm*sqrt((n-m)*(n-m-1))
                #px+= 0.5*f*(C(n-1,m+1)*Vnm+S(n-1,m+1)*Wnm)*sqrt((n-m)*(n-m-1))
                #py+= 0.5*f*(-C(n-1,m+1)*Wnm+S(n-1,m+1)*Vnm)*sqrt((n-m)*(n-m-1))
			end
            if (n>=2 && m>=2)
				cofsx[i(n-1,m-1)]+= -0.5*Vnm*sqrt((n+m)*(n+m-1))
				cofsx[i(n-1,m-1)+L]+= -0.5*Wnm*sqrt((n+m)*(n+m-1))
				cofsy[i(n-1,m-1)]+= -0.5*Wnm*sqrt((n+m)*(n+m-1))
				cofsy[i(n-1,m-1)+L]+= 0.5*Vnm*sqrt((n+m)*(n+m-1))
                #px+= 0.5*(-C(n-1,m-1)*Vnm-S(n-1,m-1)*Wnm)*sqrt((n+m)*(n+m-1))
                #py+= 0.5*(-C(n-1,m-1)*Wnm+S(n-1,m-1)*Vnm)*sqrt((n+m)*(n+m-1))
			end
            if (m==1 && n>=2)
				cofsx[i(n-1,0)]+= -Vnm*sqrt((n+1)*(n)/2)
				cofsy[i(n-1,0)]+= -Wnm*sqrt((n+1)*(n)/2)
                #px+= -C(n-1,0)*Vnm*sqrt((n+1)*(n)/2)
                #py+= -C(n-1,0)*Wnm*sqrt((n+1)*(n)/2)
			end
            if (n>=2 && n>m)
				cofsz[i(n-1,m)]+= -Vnm*sqrt((n-m)*(m+n))
				cofsz[i(n-1,m)+L]+= -Wnm*sqrt((n-m)*(m+n))
                #pz+= (-C(n-1,m)*Vnm-S(n-1,m)*Wnm)*sqrt((n-m)*(m+n))
			end
		end
	end    
    return [cofsx'; cofsy'; cofsz']
end

# ╔═╡ cdf5be9a-3557-4a6c-a98f-e3d7470afc59
fieldcoeffs(1111164.8708100126,0.0,6259542.961028692,N,R)*[cofs.G; cofs.H]

# ╔═╡ 02165855-9830-4591-9cc0-fa4e59d9a6c0
"""
Coeffient to feild vectors, This times the coeffients [G,H] give the feild.
"""
function fieldcoeffs2(x,y,z,N,R)
    #step 1 get Vs and Ws
	i(n,m)= ((n+1) * n)÷2 + m
	L= i(N,N)
	cofsx= zeros(2L)
	cofsy= zeros(2L)
	cofsz= zeros(2L)
    r= sqrt(x*x+y*y+z*z)
    r2= x*x+y*y+z*z
    Vtop= R/r
    Wtop= 0.0
    Vprev= 0.0
    Wprev= 0.0
    Vnm= Vtop
    Wnm= Wtop
    a= z*R/r2
    b= R*R/r2
    c= x*R/r2
    d= y*R/r2
    
    for m in 0:N
        sqrtmult= (1/sqrt(1))*(1/sqrt(2*m+1))
        lastsqrtmult= 0
        for n in m:N
            if(n==m)
                if(m!=0)
                    if(m==1)
                        nm_mult=1
                    else
                        nm_mult= (2*m-1)*(1/sqrt(2*m-1))*(1/sqrt(2*m))
					end
                    temp= Vtop
                    Vtop= nm_mult*(c*Vtop - d*Wtop)
                    Wtop= nm_mult*(c*Wtop + d*temp) 
                    Vprev= 0.0
                    Wprev= 0.0
                    Vnm= Vtop
                    Wnm= Wtop
				end
            else
				amult= (2*n-1)*sqrtmult*a
				bmult= (n+m-1)*(n-m-1)*sqrtmult*lastsqrtmult*b
                temp= Vnm
                Vnm= amult*Vnm - bmult*Vprev
                Vprev= temp
                temp= Wnm
                Wnm= amult*Wnm - bmult*Wprev
                Wprev= temp
                lastsqrtmult= sqrtmult
                sqrtmult= (1/sqrt(n+1-m))*(1/sqrt(n+1+m))
			end
            if (m<N && n>=m+2)
                f=  0.5*(n-m-0.58) #0.5*sqrt((n-m)*(n-m-1))
                if m==0
                    f= 0.5*sqrt(2)*sqrt((n-m)*(n-m-1))
				end
				cofsx[i(n-1,m+1)]+= f*Vnm
				cofsx[i(n-1,m+1)+L]+= f*Wnm
				cofsy[i(n-1,m+1)]+= -f*Wnm
				cofsy[i(n-1,m+1)+L]+= f*Vnm
                #px+= 0.5*f*(C(n-1,m+1)*Vnm+S(n-1,m+1)*Wnm)*sqrt((n-m)*(n-m-1))
                #py+= 0.5*f*(-C(n-1,m+1)*Wnm+S(n-1,m+1)*Vnm)*sqrt((n-m)*(n-m-1))
			end
            if (n>=2 && m>=2)
				f= 0.5*(n+m-0.5295)#0.5*sqrt((n+m)*(n+m-1))
				cofsx[i(n-1,m-1)]+= -f*Vnm
				cofsx[i(n-1,m-1)+L]+= -f*Wnm
				cofsy[i(n-1,m-1)]+= -f*Wnm
				cofsy[i(n-1,m-1)+L]+= f*Vnm
                #px+= 0.5*(-C(n-1,m-1)*Vnm-S(n-1,m-1)*Wnm)*sqrt((n+m)*(n+m-1))
                #py+= 0.5*(-C(n-1,m-1)*Wnm+S(n-1,m-1)*Vnm)*sqrt((n+m)*(n+m-1))
			end
            if (m==1 && n>=2)
				f= sqrt((n+1)*(n)/2)
				cofsx[i(n-1,0)]+= -f*Vnm
				cofsy[i(n-1,0)]+= -f*Wnm
                #px+= -C(n-1,0)*Vnm*sqrt((n+1)*(n)/2)
                #py+= -C(n-1,0)*Wnm*sqrt((n+1)*(n)/2)
			end
            if (n>=2 && n>m)
				f= (n-m)*(m+n)*lastsqrtmult
				cofsz[i(n-1,m)]+= -Vnm*f
				cofsz[i(n-1,m)+L]+= -Wnm*f
                #pz+= (-C(n-1,m)*Vnm-S(n-1,m)*Wnm)*sqrt((n-m)*(m+n))
			end
		end
	end    
    return [cofsx'; cofsy'; cofsz']
end

# ╔═╡ cebbe720-e322-47d1-a304-fc019bcb1910
begin
	npoints=10000
	rng2 = MersenneTwister(123);
	A= zeros(3npoints,2L)
	A2= zeros(3npoints,2L)
	for i in 1:npoints
		rmin= 6300_000
		hmax= 2000_000
		point= randn(rng2,3)
		point= point/norm(point)*(rmin + rand(rng2)*hmax)
		A[3i-2:3i,:]= fieldcoeffs(point...,N,R)
		A2[3i-2:3i,:]= fieldcoeffs2(point...,N,R)
	end
end

# ╔═╡ 3eec9b8c-803b-41f1-ba4c-d1c1085cda02
xm=[cofs.G; cofs.H]

# ╔═╡ 4b9ff1ba-afaf-40be-8089-59b72213e23a
bm=A*xm

# ╔═╡ 42b2c88f-e087-467b-b66a-eaeb71faf58a
x2=pinv(A2)*bm

# ╔═╡ 3087bdbc-3011-44f9-9b64-60f97a2a7c3c
b2= A2*x2

# ╔═╡ c6739c5c-e055-4ada-b83e-7a54c86d4197
maximum(abs,b2-bm)

# ╔═╡ e3423a25-f163-46fc-95fd-beb0a97412b2
function Cfp(n,m)::Int32
	i= ((n+1) * n)÷2 + m
	if i<=0 || i>length(cofs.G) || m>n
		return Int32(0)
	else
		return round(Int32,cofs.G[i]*10*256)
	end
end

# ╔═╡ 1a0cec6a-6769-4205-95ac-7c6b1fa08285
function Sfp(n,m)::Int32
	i= ((n+1) * n)÷2 + m
	if i<=0 || i>length(cofs.H) || m>n
		return Int32(0)
	else
		return round(Int32,cofs.H[i]*10*256)
	end
end

# ╔═╡ 39441dbb-b297-4501-9f26-846c0a0e348a
fpmul(x::Int32,y::Int32,bits)= Int32(widemul(x,y)>>bits)

# ╔═╡ 32cb01cf-8e13-4637-9717-27ae8bc647e0
fpmul(x::Int32,y::Int32)= Int32(widemul(x,y)>>16)

# ╔═╡ 02379a45-656c-4470-9985-759131688704
function recsqrt(i,bits)
	round(Int32,1/sqrt(i)*(1<<bits))
end

# ╔═╡ 7ed43d96-ecc9-4216-8942-fc893258f301
function recsqrt(i)
	round(Int32,1/sqrt(i)*(1<<16))
end

# ╔═╡ 840e7244-e461-4cc9-84d2-31cd0215c955
string(recsqrt.(1:(2*N+3)))

# ╔═╡ 4b7256e5-ded5-4b22-99c3-399972bc4b69
fullmagfield(1111164.8708100126,0.0,6259542.961028692,N,R,C,S)

# ╔═╡ 62755485-e904-41ba-b863-e5f043e4d5c0
rng = MersenneTwister(1234);

# ╔═╡ f7516052-cfaf-4fa7-b217-073f80c23dd7
rand(rng)

# ╔═╡ 23405a83-82fd-4452-acdd-05efd005c584
10_0000

# ╔═╡ 471fda84-fe93-4dde-87b5-138973675589
struct Model
	sqrt_lut::Vector{UInt16}
	c_lut::Vector{UInt32}
	s_lut::Vector{UInt32}
end

# ╔═╡ 3fe07f50-52ec-4a76-8fc0-6aed5fe415b2
function fullmagfieldfp(x,y,z,N,R,bits,bits2)
    #step 1 get Vs and Ws
    px::Int32=0
    py::Int32=0
    pz::Int32=0
    r= sqrt(x*x+y*y+z*z)
    r2= x*x+y*y+z*z
    Vtop::Int32= round(R/r*(1<<bits))
    Wtop::Int32=0
    Vprev::Int32=0
    Wprev::Int32=0
    Vnm::Int32= Vtop
    Wnm::Int32= Wtop
	scalea=R/r2*(1<<bits)
    a::Int32= round(z*scalea)
    b::Int32= round(R*scalea)
    c::Int32= round(x*scalea)
    d::Int32= round(y*scalea)
    
    for m::Int32 in 0:N+1
        sqrtmult::Int32= recsqrt(2*m+1,bits)
        lastsqrtmult::Int32= 0
        for n::Int32 in m:N+1
            if(n==m)
                if(m!=0)
					local nm_mult::Int32
                    if(m==1)
                        nm_mult= 1<<bits
                    else
                        nm_mult= fpmul(Int32(2*m-1)*recsqrt(2*m-1,bits), recsqrt(2*m,bits),bits)
					end
                    temp= Vtop
                    Vtop= fpmul(nm_mult,Int32((widemul(c,Vtop) - widemul(d,Wtop))>>bits),bits)
                    Wtop= fpmul(nm_mult,Int32((widemul(c,Wtop) + widemul(d,temp))>>bits),bits)
                    Vprev= 0
                    Wprev= 0
                    Vnm= Vtop
                    Wnm= Wtop
				end
            else
				amult::Int32= fpmul(sqrtmult*Int32(2*n-1),a,bits)
				bmult::Int32= fpmul(fpmul(Int32(n+m-1)*sqrtmult,Int32(n-m-1)*lastsqrtmult,bits), b,bits)
                temp= Vnm
                Vnm= (widemul(amult,Vnm) - widemul(bmult,Vprev))>>bits
                Vprev= temp
                temp= Wnm
                Wnm= (widemul(amult,Wnm) - widemul(bmult,Wprev))>>bits
                Wprev= temp
                lastsqrtmult= sqrtmult
				sqrtmult= fpmul(recsqrt(n+1-m,bits),recsqrt(n+1+m,bits),bits)
			end
			local f::Int32
            if (m<N && n>=m+2)
                if m==0
                    f= fpmul(Int32(n)*recsqrt(n,bits), Int32(n-1)*recsqrt((n-1)*2,bits),bits)
				else
					f= widemul(Int32(n-m)*recsqrt(n-m,bits), Int32(n-m-1)*recsqrt(n-m-1,bits))>>(bits+1)
				end
                px+= widemul(f,Int32((widemul(Cfp(n-1,m+1),Vnm)+widemul(Sfp(n-1,m+1),Wnm))>>bits2))>>bits
				py+= widemul(f,Int32((widemul(-Cfp(n-1,m+1),Wnm)+widemul(Sfp(n-1,m+1),Vnm))>>bits2))>>bits
			end
            if (n>=2 && m>=2)
				f= widemul(Int32(n+m)*recsqrt(n+m,bits), Int32(n+m-1)*recsqrt(n+m-1,bits))>>(bits+1)
				px+= widemul(f,Int32((widemul(-Cfp(n-1,m-1),Vnm)+widemul(-Sfp(n-1,m-1),Wnm))>>bits2))>>bits
                py+= widemul(f,Int32((widemul(-Cfp(n-1,m-1),Wnm)+widemul(Sfp(n-1,m-1),Vnm))>>bits2))>>bits
			end
            if (m==1 && n>=2)
				f= widemul(Int32(n+1)*recsqrt(n+1,bits), Int32(n)*recsqrt(n*2,bits))>>bits
				f= fpmul(f,-Cfp(n-1,0),bits2)
                px+= fpmul(f,Vnm,bits)
                py+= fpmul(f,Wnm,bits)
			end
            if (n>=2 && n>m)
				f= widemul(Int32(n-m)*recsqrt(n-m,bits), Int32(n+m)*recsqrt(n+m,bits))>>bits
				pz+= widemul(f,Int32((widemul(-Cfp(n-1,m),Vnm)+widemul(-Sfp(n-1,m),Wnm))>>bits2))>>bits
			end
		end
	end    
    return ((px,py,pz) .>> (bits-bits2)) ./ 10
end

# ╔═╡ 39ecdfe9-4c06-4d8d-8368-db710d3ba13f
function fullmagfieldfp2(x,y,z,N,R,bits,bits2)
    #step 1 get Vs and Ws
    px::Int32=0
    py::Int32=0
    pz::Int32=0
    r= sqrt(x*x+y*y+z*z)
    r2= x*x+y*y+z*z
    Vtop::Int32= round(R/r*(1<<bits))
    Wtop::Int32=0
    Vprev::Int32=0
    Wprev::Int32=0
    Vnm::Int32= Vtop
    Wnm::Int32= Wtop
	scalea=R/r2*(1<<bits)
    a::Int32= round(z*scalea)
    b::Int32= round(R*scalea)
    c::Int32= round(x*scalea)
    d::Int32= round(y*scalea)
    
    for m::Int32 in 0:N+1
        sqrtmult::Int32= recsqrt(2*m+1)
        lastsqrtmult::Int32= 0
        for n::Int32 in m:N+1
            if(n==m)
                if(m!=0)
					local nm_mult::Int32
                    if(m==1)
                        nm_mult= 1<<bits
                    else
                        nm_mult= fpmul(Int32(2*m-1)*recsqrt(2*m-1), recsqrt(2*m))
					end
                    temp= Vtop
                    Vtop= fpmul(nm_mult,fpmul(c,Vtop) - fpmul(d,Wtop))
                    Wtop= fpmul(nm_mult,fpmul(c,Wtop) + fpmul(d,temp))
                    Vprev= 0
                    Wprev= 0
                    Vnm= Vtop
                    Wnm= Wtop
				end
            else
				amult::Int32= fpmul(sqrtmult*Int32(2*n-1),a)
				bmult::Int32= fpmul(fpmul(Int32(n+m-1)*sqrtmult,Int32(n-m-1)*lastsqrtmult), b)
                temp= Vnm
                Vnm= fpmul(amult,Vnm) - fpmul(bmult,Vprev)
                Vprev= temp
                temp= Wnm
                Wnm= fpmul(amult,Wnm) - fpmul(bmult,Wprev)
                Wprev= temp
                lastsqrtmult= sqrtmult
				sqrtmult= fpmul(recsqrt(n+1-m),recsqrt(n+1+m))
			end
			local f::Int32
            if (m<N && n>=m+2)
                if m==0
                    f= fpmul(Int32(n)*recsqrt(n), Int32(n-1)*recsqrt((n-1)*2))
				else
					f= fpmul(Int32(n-m)*recsqrt(n-m), Int32(n-m-1)*recsqrt(n-m-1))>>1
				end
                px+= fpmul(f,fpmul(Cfp(n-1,m+1),Vnm)+fpmul(Sfp(n-1,m+1),Wnm) )
				py+= fpmul(f,fpmul(-Cfp(n-1,m+1),Wnm)+fpmul(Sfp(n-1,m+1),Vnm) )
			end
            if (n>=2 && m>=2)
				f= fpmul(Int32(n+m)*recsqrt(n+m), Int32(n+m-1)*recsqrt(n+m-1))>>1
				px+= fpmul(f,fpmul(-Cfp(n-1,m-1),Vnm)+fpmul(-Sfp(n-1,m-1),Wnm) )
                py+= fpmul(f,fpmul(-Cfp(n-1,m-1),Wnm)+fpmul(Sfp(n-1,m-1),Vnm) )
			end
            if (m==1 && n>=2)
				f= fpmul(Int32(n+1)*recsqrt(n+1), Int32(n)*recsqrt(n*2))
				f= fpmul(f,-Cfp(n-1,0))
                px+= fpmul(f,Vnm)
                py+= fpmul(f,Wnm)
			end
            if (n>=2 && n>m)
				f= fpmul(Int32(n-m)*recsqrt(n-m), Int32(n+m)*recsqrt(n+m))
				pz+= fpmul(f,fpmul(-Cfp(n-1,m),Vnm)+fpmul(-Sfp(n-1,m),Wnm) )
			end
		end
	end    
    return ((px,py,pz) .>> (8)) ./ 10
end

# ╔═╡ 8e84e1f1-8ae4-4fda-ad8f-ef1b62d0680d
fullmagfieldfp2(1111164.8708100126,0.0,6259542.961028692,N,R,16,8)

# ╔═╡ 4b5e70f0-3db8-4202-8a11-67bae83ae023
fullmagfieldfp2(1128529.6885767058,0.0,6358023.736329913,N,R,16,8)

# ╔═╡ 3c7bc343-bb3a-43e7-86ca-0178a9c52c56
for i in 1:100000
	rmin= 6300_000
	hmax= 2000_000
	point= randn(rng,3)
	point= point/norm(point)*(rmin + rand(rng)*hmax)
	floatb= fullmagfield(point...,N,R)
	fixb= fullmagfieldfp2(point...,N,R,16,8)
	@assert norm(floatb .- fixb)<7
end

# ╔═╡ f2455a84-abcb-4d27-95e3-ff5cac1f3d04
function fullmagfieldfp3(x,y,z,N,R,bits,bits2)
    #step 1 get Vs and Ws
    px::Int32=0
    py::Int32=0
    pz::Int32=0
    r= sqrt(x*x+y*y+z*z)
    r2= x*x+y*y+z*z
    Vtop::Int32= round(R/r*(1<<bits))
    Wtop::Int32=0
    Vprev::Int32=0
    Wprev::Int32=0
    Vnm::Int32= Vtop
    Wnm::Int32= Wtop
	scalea=R/r2*(1<<bits)
    a::Int32= round(z*scalea)
    b::Int32= round(R*scalea)
    c::Int32= round(x*scalea)
    d::Int32= round(y*scalea)
    
    for m::Int32 in 0:N+1
        sqrtmult::Int32= recsqrt(2*m+1)
        lastsqrtmult::Int32= 0
        for n::Int32 in m:N+1
            if(n==m)
                if(m!=0)
					local nm_mult::Int32
                    if(m==1)
                        nm_mult= 1<<bits
                    else
                        nm_mult= fpmul(Int32(2*m-1)*recsqrt(2*m-1), recsqrt(2*m))
					end
                    temp= Vtop
                    Vtop= fpmul(nm_mult,fpmul(c,Vtop) - fpmul(d,Wtop))
                    Wtop= fpmul(nm_mult,fpmul(c,Wtop) + fpmul(d,temp))
                    Vprev= 0
                    Wprev= 0
                    Vnm= Vtop
                    Wnm= Wtop
				end
            else
				amult::Int32= fpmul(sqrtmult*Int32(2*n-1),a)
				bmult::Int32= fpmul(fpmul(Int32(n+m-1)*sqrtmult,Int32(n-m-1)*lastsqrtmult), b)
                temp= Vnm
                Vnm= fpmul(amult,Vnm) - fpmul(bmult,Vprev)
                Vprev= temp
                temp= Wnm
                Wnm= fpmul(amult,Wnm) - fpmul(bmult,Wprev)
                Wprev= temp
                lastsqrtmult= sqrtmult
				sqrtmult= fpmul(recsqrt(n+1-m),recsqrt(n+1+m))
			end
			local f::Int32
            if (m<N && n>=m+2)
                if m==0
                    f= fpmul(Int32(n)*recsqrt(n), Int32(n-1)*recsqrt((n-1)*2))
				else
					f= fpmul(Int32(n-m)*recsqrt(n-m), Int32(n-m-1)*recsqrt(n-m-1))>>1
				end
                px+= fpmul(f,fpmul(Cfp(n-1,m+1),Vnm)+fpmul(Sfp(n-1,m+1),Wnm) )
				py+= fpmul(f,fpmul(-Cfp(n-1,m+1),Wnm)+fpmul(Sfp(n-1,m+1),Vnm) )
			end
            if (n>=2 && m>=2)
				f= fpmul(Int32(n+m)*recsqrt(n+m), Int32(n+m-1)*recsqrt(n+m-1))>>1
				px+= fpmul(f,fpmul(-Cfp(n-1,m-1),Vnm)+fpmul(-Sfp(n-1,m-1),Wnm) )
                py+= fpmul(f,fpmul(-Cfp(n-1,m-1),Wnm)+fpmul(Sfp(n-1,m-1),Vnm) )
			end
            if (m==1 && n>=2)
				f= fpmul(Int32(n+1)*recsqrt(n+1), Int32(n)*recsqrt(n*2))
				f= fpmul(f,-Cfp(n-1,0))
                px+= fpmul(f,Vnm)
                py+= fpmul(f,Wnm)
			end
            if (n>=2 && n>m)
				f= fpmul(Int32(n-m)*recsqrt(n-m), Int32(n+m)*recsqrt(n+m))
				pz+= fpmul(f,fpmul(-Cfp(n-1,m),Vnm)+fpmul(-Sfp(n-1,m),Wnm) )
			end
		end
	end    
    return ((px,py,pz) .>> (8)) ./ 10
end

# ╔═╡ 78feafbd-f11b-48b1-b145-a14dc1fc1b3a
recsqrt(1)

# ╔═╡ Cell order:
# ╠═1764fb52-c62d-11eb-0b88-9f632da0c4f7
# ╠═57f856ab-c04e-46a6-883d-44eb39c248fe
# ╠═d2647b9f-5c7d-4f7c-a559-204ad8081ad3
# ╠═22dedee1-7710-41de-82d8-453ffd9c798a
# ╠═2ec9b6be-2b88-4360-9585-893e01a8aa52
# ╠═0a0477b3-ff5f-4c13-9168-57810ee2549f
# ╠═85e7c3fc-0dd8-4b9d-a224-f6a97db2f210
# ╠═5eb05e2f-08f6-49c7-aede-4bc8a4f885e4
# ╠═db18d18c-1937-457c-a18a-c55a5f6667db
# ╠═840e7244-e461-4cc9-84d2-31cd0215c955
# ╠═3c8e550d-e407-499c-a48a-d4a31351c1be
# ╠═f60ba0ac-aab6-4fe8-a281-d50440cf09bf
# ╠═8c5baa99-7043-4658-b1bf-c4eddde7a734
# ╠═b0a5031b-6890-42eb-9655-52d57eddc1c0
# ╠═16f5d989-edf6-4c18-aa95-963ed934a33a
# ╟─63cf4e77-cda9-442f-b78d-f8feadadfe81
# ╠═bfa7f0f0-60c1-49de-b337-252548ee34e3
# ╠═048c58bc-b114-40e1-9506-cac194941cc1
# ╠═1320035e-4ffd-45c1-a502-871fc03bbf19
# ╠═cdf5be9a-3557-4a6c-a98f-e3d7470afc59
# ╠═e24a3ce2-93d6-4e2b-85fb-008edc98ab31
# ╠═4c60844a-bba9-4192-870f-a2a6c80c6f99
# ╠═02165855-9830-4591-9cc0-fa4e59d9a6c0
# ╠═c6739c5c-e055-4ada-b83e-7a54c86d4197
# ╠═cebbe720-e322-47d1-a304-fc019bcb1910
# ╠═3eec9b8c-803b-41f1-ba4c-d1c1085cda02
# ╠═4b9ff1ba-afaf-40be-8089-59b72213e23a
# ╠═42b2c88f-e087-467b-b66a-eaeb71faf58a
# ╠═3087bdbc-3011-44f9-9b64-60f97a2a7c3c
# ╠═e3423a25-f163-46fc-95fd-beb0a97412b2
# ╠═1a0cec6a-6769-4205-95ac-7c6b1fa08285
# ╠═39441dbb-b297-4501-9f26-846c0a0e348a
# ╠═32cb01cf-8e13-4637-9717-27ae8bc647e0
# ╠═02379a45-656c-4470-9985-759131688704
# ╠═7ed43d96-ecc9-4216-8942-fc893258f301
# ╠═4b7256e5-ded5-4b22-99c3-399972bc4b69
# ╠═8e84e1f1-8ae4-4fda-ad8f-ef1b62d0680d
# ╠═4b5e70f0-3db8-4202-8a11-67bae83ae023
# ╠═62755485-e904-41ba-b863-e5f043e4d5c0
# ╠═f7516052-cfaf-4fa7-b217-073f80c23dd7
# ╠═23405a83-82fd-4452-acdd-05efd005c584
# ╠═3c7bc343-bb3a-43e7-86ca-0178a9c52c56
# ╟─471fda84-fe93-4dde-87b5-138973675589
# ╟─3fe07f50-52ec-4a76-8fc0-6aed5fe415b2
# ╟─39ecdfe9-4c06-4d8d-8368-db710d3ba13f
# ╟─f2455a84-abcb-4d27-95e3-ff5cac1f3d04
# ╠═78feafbd-f11b-48b1-b145-a14dc1fc1b3a
# ╠═8fdf5004-b7c5-4575-b455-2ec40748f355
# ╠═927a8e03-2db7-4d35-b391-5d5cd13f6ea4
