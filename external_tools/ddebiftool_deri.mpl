ddebiftool_deri:=module()
	description "directional derivatives of f(xx,p) with respect to xx",
	"use sd_dirderi for state-dependent delays and dirderi for constant delays";
	uses CodeGeneration, LinearAlgebra, VectorCalculus;
	export dirderi, sd_dirderi;
	option package;
	sd_dirderi:=proc(f0::Vector, maxorder::posint)
		local dim,xsubs,f0sub,hf,f,dsubs,i,i1,fcndefs,cstr,fcnstr;
		dim:=Dimension(f0); # dimension of system;
		# replace x by equilibrium+h*deviation
		xsubs:={seq(x[jj]=x0[jj]+hf*v[jj],jj=1..dim)} ;
		f0sub:=eval(f0,xsubs);
		hf:=t->h;
		f:=eval(f0sub);
		# differentiate f i times (up to maxorder)
		for i from 1 to maxorder do
  			df0||i:=eval(diff(f,seq(h,jj=1..i)),h=0);
		od;
		# create substitutes of derivatives and indices using name concatenations
		dsubs:={seq(v[jj]=v0||jj,jj=1..dim),seq(x0[jj]=x||jj,jj=1..dim)};
		for i1 from 1 to maxorder do
	   		dsubs:=dsubs union {seq((D@@i1)(x0[jj])=0,jj=1..dim)};
   			dsubs:=dsubs union {seq((D@@i1)(v[jj])=v||i1||jj,jj=1..dim)};
		od;
		# revert xj back to x[j]
		for i1 from 1 to dim do
 			x||i1:=unapply(x[i1],t);
		od;
		# insert name concatenations into derivatives of f
		for i1 from 1 to maxorder do
 			df||i1:=eval(df0||i1,dsubs);
		od;
		# create header of matlab file (with matlab function handles vkj)
		fcndefs:=sprintf("function y=dirderi(i,x,p,vinp)\nsel=@(x0,y0)bsxfun(@(x,y)x(y),x0,y0);\n");
		for i from 0 to maxorder do
			for i1 from 1 to dim do
  				s||i||i1:=sprintf("v%d%d=@(t)sel(vinp(%d,t),%d);\n",i,i1,i,i1);
 				fcndefs:=cat(fcndefs,s||i||i1);
 			od;
		od;
		# create matlab code for derivatives using CodeGeneration
		for i from 1 to maxorder do
	 		cm||i:=Matlab(df||i,optimize,resultname=y,output=string);
		od;
		# wrap derivatives into a switch statement
		cstr:=sprintf("switch i\n");
		for i from 1 to maxorder do
 			cstr:=cat(cstr,sprintf("case %d\n",i),cm||i);
		od;
		cstr:=cat(cstr,sprintf("end\nend\n"));
		fcnstr:=cat(fcndefs,cstr);
		return fcnstr, [seq(df||jj,jj=1..maxorder)];
	end proc;
	dirderi:=proc(f0::Vector, ndelays::posint,maxorder::posint)
		local dim,xsubs,f,i,fcndefs,cstr,fcnstr;
		dim:=Dimension(f0); # dimension of system
		# replace x by equilibrium+h*deviation
		xsubs:={seq(seq(x[jj,kk]=x[jj]+h*v[jj,kk],jj=1..dim),kk=1..ndelays+1)};
		f:=eval(f0,xsubs);
		# differentiate f i times (up to maxorder)
		for i from 1 to maxorder do
	  		df||i:=eval(diff(f,seq(h,jj=1..i)),h=0);
		od;
		# create header of matlab file (with matlab function handles vkj)
		fcndefs:=sprintf("function y=dirderi(i,x,v,p)\n");
		# create matlab code for derivatives using CodeGeneration
		for i from 1 to maxorder do
	 		cm||i:=Matlab(df||i,optimize,resultname=y,output=string);
		od;
		# wrap derivatives into a switch statement
		cstr:=sprintf("switch i\n");
		for i from 1 to maxorder do
 			cstr:=cat(cstr,sprintf("case %d\n",i),cm||i);
		od;
		cstr:=cat(cstr,sprintf("end\nend\n"));
		fcnstr:=cat(fcndefs,cstr);
		return fcnstr,[seq(df||jj,jj=1..maxorder)];
	end proc;
end module;
