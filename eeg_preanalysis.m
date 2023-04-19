%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    PREANALIZA_BS POR LOTES DE SUJETOS   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Configuration files:

% Crear una matriz  (.mat) en la que la primera columna sea el número de
% condición (p.e., 1, 2, 3) y la segunda los pulsos correspondientes (p.e.,
% 1, 128, 255), y llamarla "equicondpulso"; guardarla en el directorio raíz
% de preanálisis (por debajo de las carpetas de los sujetos preanalizados;
% p.e., c:\**\**\**\preana). Servirá para todos los sujs del
% experimento. El script abre esta matriz en el 3er bloque.

% En esa misma carpeta, copiar "electrovecinos".

% Demarcar época con margen (e.g. 50ms) al principio y al final; en
% script para fabricar ANA y GPs se va a corregir el retraso en la
% aparición del estímulo y se va a definir la ventana final teniendo esto
% en cuenta.

%  Fátima Alvarez (UAM) 04/2023: Script updated from previous' Luis Carretié (UAM)script 

% Set fieltrip last version:
% 
% restoredefaultpath 
% % addpath C:\Psicofis\fieldtrip-20221212 % añadir versión más reciente
% addpath C:\Psicofis\fieldtrip-20230128; %add with subfolders
% ft_defaults
% addpath('C:\**\**\**') %añadir carpeta con script de análisis

%% Para que funcione bien fieltrip hay que instalar algunas toolbox, 
% para ello hay que ir a Home> environment > add-ons. Instalar:

% (x) Image Processing Toolbox (images)
% (x) Optimization Toolbox (optim)
% (x) Signal Processing Toolbox (signal)
% (x) Statistics and Machine Learning Toolbox (stats)

%% 
clear all
close all
clc

%=================================================
%%   Parámetros a modificar en cada experimento   %
%=================================================

% define ruta

ruta_comun='W:\**\**\**\'; % ruta datos
direc_preana='Preanalisis\'; % ruta script preanálisis

%  define sujetos:

% sujetos={ 's13','s461509','s471511','s481512', 's491515', 's501516', 's511611'... 
%     's521612', 's531615', 's541616', 's551709', 's561711', 's571712'...
%     's581715','s591716', 's601809', 's611811', 's621812', 's631815',... 
%     's641816', 's652109', 's662111', 's672112', 's682115', 's692209', ...
%     's702211', 's712212', 's722215', 's732216', 's742309', 's752311', ...
%     's762312', 's772315', 's782316', 's802412', 's812415', ...
%     's822416', 's832509', 's842511', 's852512', 's862515', 's872516', ...
%     's882809', 's892811', 's902812', 's912815', 's922816', 's932909', ...
%     's942911', 's952912', 's962915', 's972916', 's983009', 's993011', ...
%     's1003012', 's1013015', 's1023016', 's1030109', 's1040111', 's1050112',...
%     's1060115', 's1070116', 's1080209', 's1090211', 's1100212',...
%     's421409', 's431411', 's441412'}; 

%% define filtros: 

%Filtros (ojo, se pasan 2 veces, antes y después de ICA; si se cambia de
%tipo de filtro, hacerlo las dos veces)

f1=0.01; % Paso alto
f2=30; % Paso bajo (comentar el que corresponda si no es paso de banda)

%Línea base para sustracción LB
lbini=-0.1; %En segundos
lbfin=0;

%Resto de la época (intervalo tras la aparición del estímulo)
epocapost0=0.2; %En segundos (habría que subirlo a 0.208)

% Intervalo crítico (ventana de interés -vdi-) EN PUNTOS a la hora de detectar
% interferencias (normalmente, intervalo que se someterá a análisis); poner
% 1 y length(tiempo) si vdi es toda la época preanalizada. Ver vector
% tiempo (aparece tras ejecutar un sujeto) para equivalencia ms-puntos.
inivdi=1;
finvdi=51;

%% Lectura condiciones
rutaequicondpulso=[ruta_comun,direc_preana];

cd(rutaequicondpulso)
load 'equicondpulso'
load 'electrovecinos'
ncond=length(equicondpulso(:,1));
pulso_code=equicondpulso(:,2); % Pulsos correspondientes a cada condición

% NOTA: los pulsos acaban codificándose +512 (valor 0) del valor incial que se le da
% en el script de settriggers


% ************************************************************************
%% *                           COMIENZA BUCLE                             *
% ************************************************************************

 for sujeto=1:length(sujetos) %% para debuggear quita el bucle
    close all %solo cuando debuggeo
%     sujeto= 11; % <--------------------------------------- +++++++++++++++
    
    ncond=length(equicondpulso(:,1));% es necsario repetir porque hay un ncond=ncond-1
    % dentro del bucle (sólo para Emotext), borrar en otros experimentos
    
    % Rutas
    sujeto_eeg=sujetos{sujeto};
    direc_eeg=[sujeto_eeg,'\'];
    ruta_eeg=[ruta_comun, direc_eeg];
    cd(ruta_eeg)
    fichconducta=['Output_', sujeto_eeg];
    load (fichconducta)
    
    %vector con el orden de condiciones en este sujeto
    ordencondic=e(:,3);
    clear condic
    
    % Se crea subcarpeta dentro de \preana (! crear \preana antes de la 1ª
    % ejecución)
    ruta_preana=[ruta_comun,direc_preana, sujeto_eeg]; %ruta preanálisis pasa a incluir la subcarpeta del sujeto
    mkdir (ruta_preana)
    
    % Definición épocas en ms
    epinims=lbini*1000;
    epfinms=epocapost0*1000;
    
    %% =================
    % Lectura datos EEG
    % =================

    eegfile= ls(sprintf('txt_%s.bdf',sujetos{sujeto} (2:end))); 

    % Selección y definición de canales
    cfg=[];
    cfg.dataset=eegfile;
    cfg.channel={'all' '-EXG6' '-EXG7' '-EXG8' '-Status'};
    cfg.reref = 'yes';
    cfg.refchannel = 'EXG5';
    
    % Filtro offline
    cfg.bpfilter='yes'; % Primer filtrado (sin él, ICA no funciona bien)
    cfg.bpfreq=[f1 f2];
    cfg.bpfilttype='but';
    cfg.bpfiltord=2;
    data=ft_preprocessing(cfg);
    
    % Más definición de canales
    ncaneeg=length(data.cfg.channel)-5;
    data.label{65}='vEOG1';
    data.label{66}='vEOG2';
    data.label{67}='hEOG1';
    data.label{68}='hEOG2';
    data.label{69}='nose_ref';
    
    % Leer pulsos
    event= ft_read_event(eegfile);
    event = event(find(strcmp({event.type},'STATUS')));
        
    
%     %% Cuando hay que más de dos ficheros bdf para un sujeto:
%     
%     % Primero leer los datos de eeg en dos ficheros separados (ej: data y
%     % datab4, event y eventb4)
%     eegfileb4= 'txt_792411_4block.bdf' %leer fichero a pegar
%         cfg.dataset=eegfileb4; %*************************************cambiar esto despues!
%         datab4=ft_preprocessing(cfg); 
%     cfg.dataset=eegfile;
%         
%     % Crear una variable nueva con la información de ambos ficheros
%     datos=data; %copiamos primero las cosas que no hay que modificar    
%     datos.hdr.nSamples = data.hdr.nSamples + datab4.hdr.nSamples;
%     datos.time{1,1}= [data.time{1,1}  (datab4.time{1,1}+ data.time{1,1}(end))];
%     datos.trial{1,1}= [data.trial{1,1}  datab4.trial{1,1}];
%     datos.sampleinfo(2)= data.sampleinfo(2)+ datab4.sampleinfo(2);
%     datos.cfg.trl(2)= data.cfg.trl(2)+ datab4.cfg.trl(2);
%     
%     data=datos; clear datos
%     
%     % leer datos de los triggers:
%     
%     eventb4= ft_read_event(eegfileb4);
%     for e=1:length(eventb4)
%         eventb4(e).sample= eventb4(e).sample+event(end).sample
%     end
%         
%     eventb4.sample = eventb4(:).sample + event(end).sample 
%     %copiar sample and value de eventb4 en event. tambien hay que copiar el type 'STATUS' 
% %     event(602).type= 'STATUS';
%        
%     clear eventb4
       
    %%
    % ======================================================
    % Comprobación automática de si los pulsos son correctos
    % ======================================================
    
    % Asignación de pulso teórico a cada ensayo
    pulsoteo=[];
    for i=1:length(ordencondic)
        for j=1:length(equicondpulso(:,1))
            if ordencondic(i)==equicondpulso(j,1)
                pulsoteo(i)=equicondpulso(j,2);
            end
        end
    end
    clear i j 
     
    % Detección discrepancia pulsos teóricos (trig) versus reales (event.value)
    if length(event) ~= length(pulsoteo)
        disp('****** DISCREPANCIA PULSOS TEÓRICOS - REALES ******; CORRIGIENDO....');
        for evento=1:length(event)
            vectorevent(evento)=event(evento).value; %Conversión del campo event.value a vector para facilitar siguientes pasos
        end
        if length(event)>length(pulsoteo)
            [elementodiscrepante,campodiscrepante]=setdiff(vectorevent,pulsoteo); %Detección de elemento discrepante
            event(campodiscrepante)=[]; %Eliminación de elemento discrepante
            for evento=1:length(event) %Comparación pulso teórico vs real
                dif(evento)=pulsoteo(evento)-event(evento).value;
            end
           if max(dif) ~= min(dif) % Si tras las correcciones sigue habiendo discrepancia
                plot (1:length(event),dif)
                error=msgbox ('Alguna correspondencia pulso teórico-pulso real no es correcta y no ha podido corregirse');
                return
            else
                disp('>> DISCREPANCIA CORREGIDA!! :)');
            end
        end
        
        if length(event)<length(pulsoteo)
            %Primer ensayo realmente registrado (cuando hay menos, suele ser porque se
            %empezó a registrar más tarde y se perdieron los primeros ensayos)
            primerens=length(pulsoteo)-length(event);
            pulsoteo=pulsoteo(primerens+1:length(pulsoteo));
            for evento=1:length(event) %Comparación pulso teórico vs real desde el ensayo en el que comenzó a registrarse
                dif(evento)=pulsoteo(evento)-event(evento).value;
            end
            if max(dif) ~= min(dif) % Si tras las correcciones sigue habiendo discrepancia
                plot (1:length(event),dif)
                error=msgbox ('Alguna correspondencia pulso teórico-pulso real no es correcta y no ha podido corregirse');
                return
            else
                disp('>> DISCREPANCIA CORREGIDA!! :)');
            end
        end
        clear evento dif
    else
        disp('Correspondencia entre pulsos reales y teóricos correcta');
    end
    
    % Guardar datos y pulsos leídos
    data.cfg.event=event;
    
    %% =================================
    %                ICA
    % =================================
    
    %Se quitan todos los que no son canales EEG para el ICA
    cfg=[];
    cfg.channel={'all' '-vEOG1' '-vEOG2' '-hEOG1' '-hEOG2' '-nose_ref'};
    data_64=ft_preprocessing(cfg,data);
    
    cfg=[];
    cfg.method='runica';   %Infomax
    comp = ft_componentanalysis(cfg,data_64);
    
    preana.ERP.ica.comp=comp;
    
    sujeto_eeg %%
    
    %% DETECCIÓN MANUAL DE ICs con parpadeo: se plotean los datos y la topografía de cada componente y se aceptan/descartan manualmente
    % ================
    % % Si se quiere ACEPTAR un componente: pulsar la barra espaciadora
    % % Si se quiere ELIMINAR un componente: pulsar el botón del ratón
    
    load davg_biosemi64
    figure('Position', get(0, 'Screensize'))
    contador1=1;
    badica=[];
    vEOG=data.trial{1}(65,:)-data.trial{1}(66,:);
    hEOG=data.trial{1}(67,:)-data.trial{1}(68,:);
    for ic=1:size(comp.trial{1},1)
        
        C=comp.trial{1}(ic,:);
        Rv=corrcoef(vEOG,C);
        Rh=corrcoef(hEOG,C);
        
        % Topography
        subplot(6,2,[2 4 6])
        topo_ica=comp.topo(:,ic);
        davg2=davg;
        davg2.avg=repmat(topo_ica,[1 size(davg.avg,2)]);
        cfg=[];
        cfg.layout = 'biosemi64.lay';
        cfg.interactive='yes';
        cfg.figure='gcf'; %
        cfg.style='straight';%
        %cfg.zlim=[-10 10];
        ft_topoplotER(cfg,davg2)      % IC topography
        title(['IC #  ' num2str(ic)],'FontWeight','bold','FontSize',16)
        
        % Señal+EOG: para comparar la señal del IC con la del EOG
        subplot(6,2,[9 10])
        plot(resample(comp.time{1},1,10),resample(comp.trial{1}(ic,:),1,10))
        % Resample para plotear con menos resolución para que no tarde tanto
        title(['EEG' num2str(ic)])
        %set(gca,'XLim',[0 1000])
        set(gca,'YLim',[-50 50])
        
        subplot(6,2,[7 8])
        plot(resample(comp.time{1},1,10),resample(data.trial{1}(65,:)-data.trial{1}(66,:),1,10),'r')
        title(['EOGv; corr. con este IC = ' num2str((round(Rv(1,2)*1000))/1000)])
        %set(gca,'XLim',[0 1000])
        set(gca,'YLim',[-1000 1000])
        
        subplot(6,2,[11 12])
        plot(resample(comp.time{1},1,10),resample(data.trial{1}(67,:)-data.trial{1}(68,:),1,10),'r')
        title(['EOGh; corr. con este IC = ' num2str((round(Rh(1,2)*1000))/1000)])
        %set(gca,'XLim',[0 1000])
        set(gca,'YLim',[-1000 1000])
        
        keydown = waitforbuttonpress;
        if (keydown == 0)
            badica(contador1)=ic;  % Se guardan en una variable los ICs eliminados
            contador1=contador1+1;
        end
        
    end
    close(clf)
    clear contador
    
    
    % DETECCIÓN AUTOMÁTICA ICs de parpadeo: se marcan como badICA los ICs que correlacionan con EOGv
    % =============================================================================================
    
%     load davg_biosemi64 %% da problemas, hay que subir esta línea (borrar esto)
%     contador1=1;
%     badica=[];
%     vEOG=data.trial{1}(65,:)-data.trial{1}(66,:);
%     hEOG=data.trial{1}(67,:)-data.trial{1}(68,:);
%     for ic=1:size(comp.trial{1},1)
%         C=comp.trial{1}(ic,:);
%         % Correlaciones en el 1er, 2º y 3er tercio del componente:
%         % Es imp. que correlación sea alta en todos los tramos,
%         % y no en un punto concreto por una alta correlación en una
%         % interferencia, p.e.
%         Rv1=corrcoef(vEOG(1:round(length(vEOG)/3)),C(1:round(length(vEOG)/3)));
%         Rv2=corrcoef(vEOG(round(length(vEOG)/3)+1:round(length(vEOG)/3)*2),C(round(length(vEOG)/3)+1:round(length(vEOG)/3)*2));
%         Rv3=corrcoef(vEOG(round((length(vEOG)/3)*2)+1:length(vEOG)),C(round((length(vEOG)/3)*2)+1:length(vEOG)));
%         
%         Rv123=[Rv1(1,2),Rv2(1,2),Rv3(1,2)];
%        if mean(Rv123)>0.25 & min (Rv123)>0.2  %Ojo, modificar esto si se detecta más de un badica
%             badica(contador1)=ic;  % Se guardan en una variable los ICs eliminados
%             contador1=contador1+1;
%         end
%     end
%     close(clf)
% %     clear contador
%     clear contador1

%     % ¿comprobar badica?
%     if sum(badica) > 1 && length (badica)< 2
%         ic=badica
%         figure
%         topo_ica=comp.topo(:,ic);
%         davg2=davg;
%         davg2.avg=repmat(topo_ica,[1 size(davg.avg,2)]);
%         cfg=[];
%         cfg.layout = 'biosemi64.lay';
%         cfg.interactive='yes';
%         cfg.figure='gcf';
%         cfg.style='straight';
%         ft_topoplotER(cfg,davg2)      % IC topography
%         title(['IC #  ' num2str(ic)],'FontWeight','bold','FontSize',16)
%     else
%     end

% %% comprobar badICA no automático
%         ic=badica
%         figure
%         topo_ica=comp.topo(:,ic);
%         davg2=davg;
%         davg2.avg=repmat(topo_ica,[1 size(davg.avg,2)]);
%         cfg=[];
%         cfg.layout = 'biosemi64.lay';
%         cfg.interactive='yes';
%         cfg.figure='gcf';
%         cfg.style='straight';
%         ft_topoplotER(cfg,davg2)      % IC topography
%         title(['IC #  ' num2str(ic)],'FontWeight','bold','FontSize',16)


%% comprobar badICA no automático

        ic=badica 
        C=comp.trial{1}(ic,:);
        Rv=corrcoef(vEOG,C);
        Rh=corrcoef(hEOG,C);
        
        % Topography
        figure
        subplot(6,2,[2 4 6])
        topo_ica=comp.topo(:,ic);
        davg2=davg;
        davg2.avg=repmat(topo_ica,[1 size(davg.avg,2)]);
        cfg=[];
        cfg.layout = 'biosemi64.lay';
        cfg.interactive='yes';
        cfg.figure='gcf'; 
        cfg.style='straight';
        ft_topoplotER(cfg,davg2)      % IC topography
        title(['IC #  ' num2str(ic)],'FontWeight','bold','FontSize',16)
        
        % Señal+EOG: para comparar la señal del IC con la del EOG
        subplot(6,2,[9 10])
        plot(resample(comp.time{1},1,10),resample(comp.trial{1}(ic,:),1,10))
        title(['EEG' num2str(ic)])
        set(gca,'YLim',[-50 50])
        
        subplot(6,2,[7 8])
        plot(resample(comp.time{1},1,10),resample(data.trial{1}(65,:)-data.trial{1}(66,:),1,10),'r')
        title(['EOGv; corr. con este IC = ' num2str((round(Rv(1,2)*1000))/1000)])
        set(gca,'YLim',[-1000 1000])
        
        subplot(6,2,[11 12])
        plot(resample(comp.time{1},1,10),resample(data.trial{1}(67,:)-data.trial{1}(68,:),1,10),'r')
        title(['EOGh; corr. con este IC = ' num2str((round(Rh(1,2)*1000))/1000)])
        set(gca,'YLim',[-1000 1000])


%% ----------------------------------------- avanzar tras comprobar BADICA

    % ======================================================================
    % Se restan los componentes rechazados de los datos EEG y se aplican de
    % nuevo los filtros
    % ======================================================================
    
    cfg=[];
    cfg.component=badica;
    data_ica_prov = ft_rejectcomponent(cfg,comp,data_64);
    
    % Filtro offline
    cfg.bpfilter='yes'; % Segundo filtrado (la resta de badica vu elve a "descolocar" la señal)
    cfg.bpfreq=[f1 f2];
    cfg.bpfilttype='but';
    cfg.bpfiltord=2;
    data_ica=ft_preprocessing(cfg,data_ica_prov);
    
    data2=data;
    for i=1:length(data_64.label)
        data2.trial{1}(i,:)=data_ica.trial{1}(i,:);   % datos limpios después del ICA
    end
    
    preana.ERP.ica.badica=badica;
    preana.ERP.ica.comp=comp;
    preana.ERP.datos.prepro_ica=data2;
    
    % ===================================================
    %         Epoqueado y sustracción línea base
    % ===================================================
    ncond=ncond-1; %Sólo para Emotext: quitamos la última condición (código 250), que es la que marca el cambio de color del punto de fijación
    for cond=1:ncond
        cfg=[];
        cfg.trialfun = 'lab8_trialfun_corr'; %% tambien a error porque no esta en la ruta que tiene que estar
        cfg.prestim = lbfin-lbini;
        cfg.poststim = epocapost0;
        cfg.event = event;
        cfg.fsample = data2.fsample;
        cfg.code = pulso_code(cond);
        cfg.resp = ones(length(event),1);% Se admiten todas las épocas (no hubo conducta pero es necesario indicarlo)
        cfg2 = ft_definetrial(cfg);
        
        % Eliminar ensayo si el punto de inicio es negativo (pasa en algún
        % suj/cond)
        if cfg2.trl(1,1)<0
            cfg2.trl(1,:)=[];
        end
        
        cfg3=[];
        cfg3.trl=cfg2.trl;
        data3=ft_redefinetrial(cfg3,data2);
        
        %Sustracción línea base
        cfg=[];
        cfg.demean='yes';   %subtracts the mean of the specified baseline window
        cfg.baselinewindow = [lbini lbfin];
%         cfg.detrend ='yes'; %********************************** detrend +
        datos{cond}=ft_preprocessing(cfg,data3);
    end
    
    ntrials=[];
    for cond=1:ncond
        ntrials(cond)=length(datos{cond}.trial);
    end
    
    preana.ERP.datos.prepro_ica_epocas=datos;
    preana.ERP.ntrials=ntrials;
    
    
    % =============================================
    % Detección automática de canales a interpolar
    % =============================================
    % Funciona bien, excepto cuando el malo es Iz (can 28); en este caso
    % "arrastra" a sus vecinos (O1, Oz y O2) a interpolación. Al repasar
    % los badchannel de todos los sujetos, estar pendientes de si están
    % incluidos estos 4 canales: en ese caso será Iz el único que esté mal; 
    % probablemente. Re-preanalizar ese sujeto con script no-automático.
    
    medorig=[];
    contador1=1;contador2=1;
%     cansospcorr=[];
    cansospampl=[];badch=[];
    sumamplabsvecinos=zeros(ncaneeg,ncond);
    
    for cond=1:ncond
        cfg=[];
        datm{cond}=ft_timelockanalysis(cfg,datos{cond});
        medorig(:,:,cond)=datm{cond}.avg;
    end
    
    tiempo=datm{1}.time*1000;
    
    % load electrovecinos
    
    for can=1:ncaneeg
        % Detección discrepancias significativas con vecinos
        vecinosindex=~isnan(electrovecinos(can,:));
        vecinos=electrovecinos(can,vecinosindex);
        vecinos=vecinos(2:length(vecinos)); % Los vecinos realmente empiezan en segunda columna
        
        for vecino=1:length(vecinos)
            for cond=1:ncond
                % Mecanismo 1: baja correlación
                %                 corrprov=corrcoef(medorig(can,:,cond),medorig(vecinos(vecino),:,cond)); %Se calcula las correlaciones para la condición crítica, asumiendo que si una está ruidosa, ya es canal descartable
                %                 if abs(corrprov(1,2))<0.05
                %                     cansospcorr(contador1)=can;
                %                     contador1=contador1+1;
                %                 end
                % Mecanismo 2 (funciona mejor): amplitudes absolutas muy diferentes
                sumamplabsvecinos(can,cond)=sumamplabsvecinos(can,cond)+mean(abs(medorig(vecinos(vecino),:,cond)),2);
            end
        end
        medamplabsvecinos=sumamplabsvecinos/length(vecinos);
       
    end
    
    medamplabs=squeeze(mean(abs(medorig(:,:,:)),2));
    
    for can=1:ncaneeg
        if mean(abs(medamplabs(can,:)),2)>mean(medamplabsvecinos(can,:),2)*2 %| mean(abs(medorig(can,:,:)),2)<mean(medamplabsvecinos(can,:),2)/2
            cansospampl(contador2)=can;
            contador2=contador2+1;
        end
    end
    
    badch=cansospampl;
    
    % Interpolación canales artefactados (badch)
    % ==========================================
    
    load elec1005.mat  % plantilla con nombres y coordenadas de electrodos
    elec2=elec;
    elec.pnt=[];
    elec.label={};
    elec.label{65}='vEOG1';
    elec.label{66}='vEOG2';
    elec.label{67}='hEOG1';
    elec.label{68}='hEOG2';
    elec.label{69}='nose_ref';
    for i=1:length(datm{1}.label)-5
        elec.pnt(i,1:3) = elec2.pnt(find(strcmp(elec2.label,datm{1}.label{i})==1),:);
        elec.label{i}=datm{1}.label{i};
    end
    elec.pnt(65:69,:)=NaN;   % Externos
    
    cfg=[];
    cfg.elec=elec;
    cfg.neighbourdist=550;% Para interpolar se usan todos los canales en un círculo de 5.5cm de diámetro
    cfg.method= 'distance'; %*************************************ACTUALIZACIÓN
    neighbours = ft_prepare_neighbours(cfg, datos{1});%***********ACTUALIZACIÓN

    medinterpol=medorig;
    ch_interpol={};
    for ch1=1:length(badch)
        can=badch(ch1);
        ct=1;
        chneig=[];
        for ch2=1:length(neighbours(can).neighblabel)
            chlab=neighbours(can).neighblabel(ch2);
            chneig(ct)=find(strcmp(elec.label,chlab)==1);
            ct=ct+1;
        end
        
        ct=1;
        chneig2=[];
        for i=1:length(chneig)
            if sum(badch==chneig(i))==0
                chneig2(ct)=chneig(i);
                ct=ct+1;
            end
        end
        
        medinterpol(can,:,:)=mean(medorig([chneig2],:,:),1);
        ch_interpol{ch1}=chneig2;
    end
    
    datos2=datos;
    for cond=1:length(datos)
        for ens=1:length(datos{cond}.trial)
            for ch1=1:length(badch)
                can=badch(ch1);
                chneig2=ch_interpol{ch1};
                datos2{cond}.trial{ens}(can,:)=mean(datos{cond}.trial{ens}([chneig2],:),1);
            end
        end
    end
    
    preana.ERP.badchannels=badch;
    
    %% =====================================
    % Rechazo automático de épocas ruidosas
    % =====================================
    
    % Si se quiere ACEPTAR un ensayo: pulsar la barra espaciadora
    % Si se quiere RECHAZAR un ensayo: pulsar el botón del ratón
    
    % Clasificación de canales (para asignarles colores distintos)
    anteriores=[1:7,33:42];
    centrales=[8:19, 32, 43:56];
    posteriores=[20:31,57:64];
    
    %Se detecta el número máximo de ensayos x condición(no siempre el nº de
    %ensayos es igual para todas las condiciones)
    for cond=1:length(datos2)
        nensxcond(cond)=length(datos2{cond}.trial);
    end
    nens=max(nensxcond);
    
    % Cálculo de la diferencia entre máximo y mínimo en cada canal y ensayo
    contador1=0;
    difintracanxcond=NaN(ncond,nens,length(datos2{1}.label));
    for cond=1:ncond
        for ens=1:length(datos2{cond}.trial)
            difintracanxcond(cond,ens,:)=max(datos2{cond}.trial{ens}(:,inivdi:finvdi),[],2)-min(datos2{cond}.trial{ens}(:,inivdi:finvdi),[],2);
            contador1=contador1+1;
            difintracan(contador1,:)=difintracanxcond(cond,ens,:);%Juntamos condiciones para media y DT
        end
    end
    
    % Umbral para detectar ensayos sospechosos (basado en media y DT)
    umbral=mean(difintracan,1)+std(difintracan,0,1).*3.5;% Reducir el multiplicador final si se quiere ser más conservador; en participantes ruidosos, mejor 3; en limpios, 3.5
    
    % Nos quedamos con los 64 canales EEG
    umbral=umbral(1:ncaneeg);
    difintracanxcond=difintracanxcond(:,:,1:ncaneeg);
    
    %Detección de las épocas en las que se sobrepasa el umbral o se sobrepasa el valor absoluto de +100, -100 µV
    sobrepasa=zeros(ncond,nens);
    canresalta=zeros(ncond,nens);
    for cond=1:ncond
        for ens=1:length(datos2{cond}.trial)
            for can=1:ncaneeg
                if difintracanxcond(cond,ens,can)>umbral(can)||max(datos2{cond}.trial{ens}(can,inivdi:finvdi),[],2)>100||min(datos2{cond}.trial{ens}(can,inivdi:finvdi),[],2)<-100
                    sobrepasa(cond,ens)=1;
                    canresalta(cond,ens)=can; % Sólo relevante para opción B
                end
            end
        end
    end
    
%     OPCIÓN A: rechazo automático de ensayos detectados en bucle anterior
    for cond=1:ncond
        contador2=1;
        for ens=1:length(datos2{cond}.trial)
            if sobrepasa(cond,ens)==1 % Comentar este if si se activa el bloque de "representación ensayos"
                badtr{cond}(contador2)=ens;
                contador2=contador2+1;
            end
        end
    end
    
%         %% OPCIÓN B-Representación ensayos - activar para opción semiautomática (se revisan
%         %visualmente los ensayos ruidosos; ver comentarios en bucles anteriores)
%         figure('Position', get(0, 'Screensize'))
%         badtr={};
%         for cond=1:ncond
%             contador2=1;
%             for ens=1:length(datos2{cond}.trial)
%                 subplot(2,1,1)           % vEOG (azul) + hEOG (rojo)
%                 plot(datos2{cond}.time{ens},datos2{cond}.trial{ens}(65,:)-datos2{cond}.trial{ens}(66,:),'b'), hold on
%                 plot(datos2{cond}.time{ens},datos2{cond}.trial{ens}(67,:)-datos2{cond}.trial{ens}(68,:),'r'), hold off
%                 set(gca,'XTick',[0:0.1:epocapost0])
%                 set(gca,'XLim',[lbini epocapost0])
%                 set(gca,'YLim',[-100 100])
%                 legend('EOGv','EOGh')
%                 title(['Condicion: ' num2str(cond) ' - Ensayo: ' num2str(ens)])
%     
%                 subplot(2,1,2)
%                 plot(datos2{cond}.time{ens},datos2{cond}.trial{ens}(anteriores,:),'color',[0.7 0 0]),hold on
%                 plot(datos2{cond}.time{ens},datos2{cond}.trial{ens}(centrales,:),'k'), hold on
%                 if canresalta(cond,ens)>0
%                     plot(datos2{cond}.time{ens},datos2{cond}.trial{ens}(posteriores,:),'color',[0 0 0.7]), hold on
%                     plot(datos2{cond}.time{ens},datos2{cond}.trial{ens}(canresalta(cond,ens),:),'color',[0 1 0],'linewidth',4), hold off
%                 else
%                     plot(datos2{cond}.time{ens},datos2{cond}.trial{ens}(posteriores,:),'color',[0 0 0.7]), hold off
%                 end
%                 set(gca,'XTick',[0:0.1:epocapost0])
%                 set(gca,'XLim',[lbini epocapost0])
%                 set(gca,'YLim',[-100 100])
%                 legend ('Ant(rojo)','Cent(negro)','Post(azul)',[num2str(canresalta(cond,ens)),' (verde)'])
%     
%                 % Si se detecta ensayo sospechoso espera a decisión, si no, avanza
%                 if sobrepasa(cond,ens)==1
%                     title ('EEG; **** barra Acepta, ratón Elimina ****', 'color','r')
%                     keydown = waitforbuttonpress;
%                     if (keydown == 0)
%                         badtr{cond}(contador2)=ens;
%                         contador2=contador2+1;
%                     end
%                 else
%                     title ('EEG')
%                     pause(0.3)
%                 end
%     
%             end
%         end
%         close(clf)
%         preana.ERP.badtrials=badtr;
%         clear nensxcond can ens difintracan contador1 contador2 sobrepasa difintracanxcond difintracan
%     
%% sigue para ambas opciones:
    
    % Nos quedamos sólo con los ensayos buenos
    % ========================================
    
    for cond=1:length(datos2)
        temp=ones(1,length(datos2{cond}.trial));
        try
            temp(badtr{cond})=0;
        end
        cfg=[];
        cfg.trials=find(temp==1);
        datos{cond}=ft_redefinetrial(cfg,datos2{cond});
    end
    preana.ERP.datos.limpios=datos;
    
    ntrials_finales=[];
    for cond=1:length(datos)
        ntrials_finales(cond)=length(datos{cond}.trial);
    end
    preana.ERP.ntrialsfinal=ntrials_finales;
    
    cd(ruta_preana)
    save preana preana -v7.3 
    exportgraphics(gca,"mybadICA.jpg","Resolution",300) % guarda el badICA para que podamos revisarlo después

    
    %% ======================================================
    % Calcula meds para los datos_limpios y gráficos finales
    % ======================================================
    
    medfin=[];
    for cond=1:length(datos)
        cfg=[];
        datm{cond}=ft_timelockanalysis(cfg,datos{cond});
        medfin(:,:,cond)=datm{cond}.avg;
    end
    
    
    % Figuras con meds finales (limpios)
    % ==================================
    rgb=rand(ncond,3);
    marcasx=[epinims:100:epfinms];
    for marx=1:length(marcasx) %Bucle para etiquetar únicamente 0 y mitad eje X
        if marcasx(marx)==0 || marcasx(marx)==epfinms/2
            etiqx{marx}=num2str(marcasx(marx));
        else
            etiqx{marx}='';
        end
    end
    
    figure ('name',sujeto_eeg,'NumberTitle','off','Position', get(0, 'Screensize'))
    posic=1;
    for can=1:64
        for cond=1:ncond
            subplot(8,8,posic)
            plot(tiempo,medfin(can,:,cond),'Color',rgb(cond,:))
            hold on
        end
        axis([tiempo(1) tiempo(length(tiempo)) -30 30])
        set(gca, 'YGrid', 'off', 'XGrid', 'on', 'XTick',marcasx,'XTickLabel',etiqx)
        title([num2str(can),' ',num2str(datm{1}.label{can})])
        posic=posic+1;
    end
    
    cd(ruta_preana)
    nomgrafmed=[sujeto_eeg,'meds.fig'];
    saveas(gcf,nomgrafmed)
    saveas(gcf,'mymeds.jpg')
	close

    % Añadir ncond en las excepciones (en N40 no, por lo que se explica al
    % comenzar el bucle:
    clearvars -except direc_preana ruta_comun sujetos f1 f2 lbini lbfin inivdi finvdi epocapost0 pulso_code equicondpulso electrovecinos
    fprintf ('Sujeto preanalizado.. siguiente \n')

 end

%% FIN
close all
