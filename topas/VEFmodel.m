classdef VEFmodel
    %Defines virtual energy fluence model, https://doi.org/10.1118/1.1543152
    %Only 2 photon sources for now, 
    %additional electron source at zs can be used to model contamination (only ~1% of fluence)

    properties
        z0
        zs
        zM
        zX
        zY
        zI
        sigma0
        sigmas
        wXI 
        wYI 
        wX
        wY
        wMX
        wMY
        P0 = 0.9;
        scaling = 1;
        EnSpectrum = 6;
        FWHM
        fluenceParam
    end
    
    methods
        function obj = VEFmodel()
        end
        
        function obj = setVEFmodel(z0,zs,zM,zX,zY,zI,sigma0,sigmas,wX,wY,wMX,wMY,EnSpectrum,P0)
            obj.z0 = z0; 
            obj.zs = zs;
            obj.zM = zM;
            obj.zX = zX;
            obj.zY = zY;
            obj.zI = zI;
            obj.sigma0 = sigma0;
            obj.sigmas = sigmas;
            obj.wX = wX;
            obj.wY = wY;
            obj.wMX = wMX;
            obj.wMY = wMY;
            obj.P0 = P0;
            obj.scaling = 1;
            obj.EnSpectrum = EnSpectrum;
        end
        
        function F = fluenceModel(obj,param,x)
            %Fluenzfunktionen (Glchg. 9 und dafür auch Glchg. 5-8) in Abhängigkeit von Parametern definieren,
            %sodass man sie als input für fitting an messdaten nehmen
            %könnte --> fit the geometrical parameters P0,sigma0,sigmaS,h0,h1,h2,h3, and h4
            obj.P0 = param(1);
            obj.scaling = param(2);
            obj.sigma0 = param(3);
            obj.sigmas= param(4);
            h0 = param(5);
            h1 = param(6);
            h2 = param(7);
            h3 = param(8);
            h4 = param(9);
            
            ind = [1:length(x)];
            x0p = min((obj.wXI.*obj.zX.*(x(:,3)-obj.z0)+2*x(:,1)*obj.zI.*(obj.z0-obj.zX))./(2.*sqrt(2).*obj.zI.*(x(:,3)-obj.zX)),(obj.wXI.*obj.zM.*(x(:,3)-obj.z0)+2.*x(:,1).*obj.zI.*(obj.z0-obj.zM))./(2.*sqrt(2).*obj.zI.*(x(:,3)-obj.zM)));
            x0m = min((obj.wXI.*obj.zX.*(x(:,3)-obj.z0)-2*x(:,1)*obj.zI.*(obj.z0-obj.zX))./(2.*sqrt(2).*obj.zI.*(x(:,3)-obj.zX)),(obj.wXI.*obj.zM.*(x(:,3)-obj.z0)-2.*x(:,1).*obj.zI.*(obj.z0-obj.zM))./(2.*sqrt(2).*obj.zI.*(x(:,3)-obj.zM)));
            y0p = min((obj.wYI.*obj.zY.*(x(:,3)-obj.z0)+2*x(:,2)*obj.zI.*(obj.z0-obj.zY))./(2.*sqrt(2).*obj.zI.*(x(:,3)-obj.zY)),(obj.wYI.*obj.zM.*(x(:,3)-obj.z0)+2.*x(:,2).*obj.zI.*(obj.z0-obj.zM))./(2.*sqrt(2).*obj.zI.*(x(:,3)-obj.zM)));
            y0m = min((obj.wYI.*obj.zY.*(x(:,3)-obj.z0)-2*x(:,2)*obj.zI.*(obj.z0-obj.zY))./(2.*sqrt(2).*obj.zI.*(x(:,3)-obj.zY)),(obj.wYI.*obj.zM.*(x(:,3)-obj.z0)-2.*x(:,2).*obj.zI.*(obj.z0-obj.zM))./(2.*sqrt(2).*obj.zI.*(x(:,3)-obj.zM)));
            F0 = ((obj.zX-obj.z0).*(obj.zY-obj.z0))./(x(:,3)-obj.z0).^2 /4.*(erf(x0p./obj.sigma0)+erf(x0m./obj.sigma0)).*(erf(y0p./obj.sigma0)+erf(y0m./obj.sigma0));
            F0(isnan(F0)) = interp1(ind(~isnan(F0)),F0(~isnan(F0)),ind(isnan(F0)), 'spline') ;
            xsp = min((obj.wXI.*obj.zX*(x(:,3)-obj.zs)+2*x(:,1)*obj.zI.*(obj.zs-obj.zX))./(2.*sqrt(2).*obj.zI.*(x(:,3)-obj.zX)),(obj.wXI.*obj.zM.*(x(:,3)-obj.zs)+2.*x(:,1)*obj.zI.*(obj.zs-obj.zM))./(2.*sqrt(2).*obj.zI.*(x(:,3)-obj.zM)));
            xsm = min((obj.wXI.*obj.zX*(x(:,3)-obj.zs)-2*x(:,1)*obj.zI.*(obj.zs-obj.zX))./(2.*sqrt(2).*obj.zI.*(x(:,3)-obj.zX)),(obj.wXI.*obj.zM.*(x(:,3)-obj.zs)-2.*x(:,1)*obj.zI.*(obj.zs-obj.zM))./(2.*sqrt(2).*obj.zI.*(x(:,3)-obj.zM)));
            ysp = min((obj.wYI.*obj.zY*(x(:,3)-obj.zs)+2*x(:,2)*obj.zI.*(obj.zs-obj.zY))./(2.*sqrt(2).*obj.zI.*(x(:,3)-obj.zY)),(obj.wYI.*obj.zM.*(x(:,3)-obj.zs)+2.*x(:,2)*obj.zI.*(obj.zs-obj.zM))./(2.*sqrt(2).*obj.zI.*(x(:,3)-obj.zM)));
            ysm = min((obj.wYI.*obj.zY*(x(:,3)-obj.zs)-2*x(:,2)*obj.zI.*(obj.zs-obj.zY))./(2.*sqrt(2).*obj.zI.*(x(:,3)-obj.zY)),(obj.wYI.*obj.zM.*(x(:,3)-obj.zs)-2.*x(:,2)*obj.zI.*(obj.zs-obj.zM))./(2.*sqrt(2).*obj.zI.*(x(:,3)-obj.zM)));
            Fs = ((obj.zX-obj.zs).*(obj.zY-obj.zs))./(x(:,3)-obj.zs).^2 ./4.*(erf(xsp/obj.sigmas)+erf(xsm./obj.sigmas)).*(erf(ysp./obj.sigmas)+erf(ysm./obj.sigmas));
            Fs(isnan(Fs)) = interp1(ind(~isnan(Fs)),F0(~isnan(Fs)),ind(isnan(Fs)), 'spline') ;
            p=sqrt(x(:,1).^2+x(:,2).^2)./(x(:,3)-obj.z0);
            F_horn = 1+p.^2.*(h0+h1.*p+h2.*p.^2+h3.*p.^3+h4.*p.^4);
            F_horn(isnan(F_horn)) = 1;
            F = obj.scaling.*(obj.P0.*F0.*F_horn + (1-obj.P0).*Fs);
            
        end
        
        function obj = setDefaultProperties(obj)
            %See phase space header Siemens PRIMUS 6MV photon beam
            obj.z0 = 0; %mm 
            obj.zs = 10; %ca. 1cm behind target?
            obj.zM = 100; %ca. 10cm behind target?
            obj.zI = 1000;%mm (1m)
            obj.zX = 350;  %ca. 35cm behind target?
            obj.zY = 400; %ca. 40cm behind target?
            obj.FWHM = 1;
            obj.sigma0 = obj.FWHM/(2*sqrt(2*log(2)));
            %--> müssen abhängig von breiten sigma0/Jaw öffnungen berechnet
            %werden
            obj.wXI = 120;
            obj.wYI = 120;
            obj.sigmas = obj.sigma0;
            obj.wX = obj.wXI*(obj.zX/obj.zI);
            obj.wY = obj.wYI*(obj.zY/obj.zI);
            obj.wMX = obj.wXI*(obj.zM/obj.zI);
            obj.wMY = obj.wYI*(obj.zM/obj.zI);
        end
        
%         function writeToFile(obj,fID)
%         %Write initial energy spectrum
%         %Muss noch richtig implementiert werden, siehe auch matRad_TopasConfig Z. 1500ff
%         fprintf(fileID,'s:So/Example/BeamEnergySpectrumType       = "Continuous"\n');
%         fprintf(fileID,'dv:So/Example/BeamEnergySpectrumValues    = %i', length(obj.EnSpectrum));
%         fprintf(fileID,'%i ',[obj.EnSpectrum.points]);
%         fprintf(fileID,' MeV\n');
%         fprintf(fileID,'uv:So/Example/BeamEnergySpectrumWeights   = %i', length(obj.EnSpectrum));
%         fprintf(fileID,'%i \n',[obj.EnSpectrum.weights./sum(obj.EnSpectrum.weights)]);
%         end
%         
        function readEnergySpectrum(obj,phaseSpaceFile,fileName,baseData)
        if phaseSpaceFile
            fID = fopen(fileName);
            l=1;
            while true
                tline = fgetl(fID);
                if contains(tline, "Energy spectrum")
                     break; 
                end   %end of file
                disp(tline);
                l = l+1;
            end
            %!würde so bis zum Ende der Datei lesen! Viell. besser einfach
            %von Datei mit 1-2 header Zeilen ausgehen statt von phase space
            %header file
            C = importdata(fileName,',',l+1);
            fclose(fID);
        end
        
        if baseData
            EnSpectrum = [baseData.machine.data(:).energySpectrum];
            obj.EnSpectrum.points = [EnSpectrum(:).energy_MeVpN];
            obj.EnSpectrum.points = [EnSpectrum(:).weight];
        end
        end
        
        function obj = fitToData(obj,filename)
            %Fluenzmodelle nutzen um Parameter an echte Fluenzdaten zu
            %fitten, Frage: Welche daten sind überhaupt vorhanden, muss
            %etwas vereinfacht werden, was sind gute Startwerte, etc.
            
            %read data from file assumed to be formatted like
            %HalcyonBeamData.xlsx
            data = {};
            data.storageInfo = spreadsheetDatastore(filename, 'Sheets',[3:11]);
            data.depthDose = fillmissing(readtable(filename,'Sheet','Open Field Depth Dose','Range','A6:H307'),'constant',0) ;
            data.latprofiles{1,1} = [1.3 5 10 20 30];
            for i=1:length(data.latprofiles{1})
            data.latprofiles{i+1,1} = fillmissing(readtable(filename,'Sheet',sprintf('Open Field Profiles at %gcm',data.latprofiles{1,1}(i)),'Range','A9:H513'),'constant',0);
            end
            data.diagprofiles = fillmissing(readtable(filename,'Sheet','Diagonal Profiles','Range','A7:F882'),'constant',0);
            
            %identify field size index -> find value in table closest to
            %calculated size or let this be class varible set by user
            %set to 4x4xcm = index 3 for now
            fieldSize = 3;
            
            %collect data for fitting
            X =[]; y = [];
            %depth data
            X = cat(1,X,[zeros(length(data.depthDose.(1)),1) zeros(length(data.depthDose.(1)),1) data.depthDose.(1)]);
            y = cat(1,y,data.depthDose.(fieldSize));
            %lateral offset profiles
            for i=1:length(data.latprofiles{1})
                %x shift
                X = cat(1,X,[data.latprofiles{i+1}.(1) zeros(length(data.latprofiles{i+1}.(1)),1) ones(length(data.latprofiles{i+1}.(1)),1).*data.latprofiles{1}(i)]);
                y = cat(1,y,data.latprofiles{i+1}.(i+1).*y(find(data.depthDose.(1)==data.latprofiles{1}(i))));
                %y shift
                X = cat(1,X,[zeros(length(data.latprofiles{i+1}.(1)),1) data.latprofiles{i+1}.(1) ones(length(data.latprofiles{i+1}.(1)),1).*data.latprofiles{1}(i)]);
                y = cat(1,y,data.latprofiles{i+1}.(i+1).*y(find(data.depthDose.(1)==data.latprofiles{1}(i))));
            end
            %diagonal offset profiles -> scale to field width since data
            %only available for 28x28cm (hardcoded to 4/divide by 7 atm)
            for i=1:size(data.diagprofiles,2)-1
                X = cat(1,X,[data.diagprofiles.(1)./7 data.diagprofiles.(1)./7 ones(length(data.diagprofiles.(1)),1).*data.latprofiles{1}(i)]);
                y = cat(1,y,data.diagprofiles.(i+1).*y(find(data.depthDose.(1)==data.latprofiles{1}(i))));
            end
            %cm to mm and shift to SAD=1m
            X = X*10;
            X(:,3) = X(:,3) + 1000;
            
            %normalize data to 0-1 interval (percentage of total dose)
            y = y./max(y);
            
            %fit fluence model
            %initial parameter guess
            param0(1) = obj.P0;
            param0(2) = obj.scaling;
            param0(3) = obj.sigma0;
            param0(4) = obj.sigmas;
            param0(5) = 0;
            param0(6) = 0;
            param0(7) = 0;
            param0(8) = 0;
            param0(9) =0;
            
            %Set constraints
            lb = [0,-Inf,0,0,-Inf,-Inf,-Inf,-Inf,-Inf]; %maybe all of them are positive? lb = [0,0,0,0,0,0,0,0];
            ub = [1,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf];
            
            options = optimoptions('lsqcurvefit','Algorithm','trust-region-reflective' ,'MaxFunctionEvaluations',1e6,'MaxIterations',1e6) %To plot convergence: 'PlotFcn','optimplotfirstorderopt')
            params = lsqcurvefit(@obj.fluenceModel,param0,X,y,lb,ub,options);
            %params2 = nlinfit(X,y,@obj.fluenceModel,param0);
            obj.fluenceParam = params;
        end
    end
end

