function find_elms_GUI

%This is a GUI in support of an algorithm for ELMs detection.
%The GUI allows the user to download data from any valid filterscope in DIII-D, without posing any constraint.
%Only fast data should be downloded as they allow more accurate estimates of the ELMs triggering time. 
%
%The algorithm is based on statistical properties and is therefore subject to false positive and missed detections.
%Its performance is determined by two parameters:
% * S-FAC, should be between 0.1 and 1. Higher values raise the threshold above which any spike is considered to be an ELM.
% * SHIFT, should be between 0 and 5 ms. Larger values discard multiple detections during an ELM [Default value is recommended]
%
%The checkbutton "autoset" automatically sets the parameter s-fac to the expected value that generally works well
%In mixed type-I/type-III elms it might be necessary to uncheck the autoset button, and set S-FAC=0.6 and SHIFT=0.2 ms
%
%Although the performance of the algorithm can be improved by tuning the parameter S-FAC, pushbuttons ADD POINT and 
%REMOVE POINT allow the user to add and remove one point, respectively, by directly clicking on the signal time trace.
%The ADD POINT pushbutton will add a point exaclty at the location indicated by the pointer, whilst the REMOVE POINT
%pushbutton will remove the closest point to the pointer.
%Multiple points can be added or removed by checking the "MULTI POINTS" checkbutton. The sequence is terminated when the last 
%added or removed point is selected with the central or the right mouse button, or by pressing anything on the keyboard.
%
%Various time windows can be analyzed by varying parameters TMIN and TMAX. When varying these two parameters data should 
%not be downloaded again for the same shot number as the GUI stores the entire dataset.
%
%Pushbutton SAVE ASCII prints out in the current directory an ascii file with the times at which ELMS are triggered. 
%The file name can be specified in the appropriate box under the SAVE ASCII pushbutton. If this is left empty, the GUI saves 
%data in a file named "shotnumber_tmin_tmax.txt".
%
%Leaving empty any of the parameters TMIN, TMAX, S-FAC, SHIFT will cause the GUI to reset its value to default.
%
%Error messages and infos on what the GUI has just done are printed in the bottom right corner of the interface. 
%When red the GUI is busy or has encountered an error, whilst when green the GUI is idle. 
%
%A. Marinoni 19/11/2014

if isempty(gcbo)
   what='incipit';
else
   what=get(gcbo,'tag');
end
strcb='find_elms_GUI';

switch what

   case 'incipit'
      
      h=figure;
      set(h,'units','normalized','position',[0.2 0.2 0.6 0.7],'name','ELMs detection v0.9','NumberTitle','off');
      set(0,'defaultuicontrolunits','normalized');
      uicontrol('position',[0.03 0.87 0.12 0.05],'style','edit','tag','shotnumber');
      uicontrol('position',[0.03 0.92 0.12 0.05],'style','text','string','shot number');
      uicontrol('position',[0.03 0.75 0.12 0.05],'style','edit','tag','filterscope','string','fs04f');
      uicontrol('position',[0.03 0.80 0.12 0.05],'style','text','string','Filterscope');
      uicontrol('position',[0.03 0.68 0.12 0.05],'callback',strcb,'tag','download','string','Download');
      uicontrol('position',[0.03 0.56 0.055 0.05],'style','edit','tag','tmin','enable','off','string','0'); 
      uicontrol('position',[0.03 0.61 0.055 0.05],'style','text','string','t min [s]','enable','off','tag','tmin_t'); 
      uicontrol('position',[0.095 0.56 0.055 0.05],'style','edit','tag','tmax','enable','off','string','4'); 
      uicontrol('position',[0.095 0.61 0.055 0.05],'style','text','string','t max [s]','enable','off','tag','tmax_t'); 
      uicontrol('position',[0.03 0.44 0.055 0.05],'style','edit','tag','sfac','string','0.2','enable','off');
      uicontrol('position',[0.03 0.49 0.055 0.05],'style','text','string','s-fac','enable','off','tag','sfac_t');
      uicontrol('position',[0.095 0.44 0.055 0.05],'style','edit','tag','shift','string','0.5','enable','off');
      uicontrol('position',[0.095 0.49 0.055 0.05],'style','text','string','shift [ms]','enable','off','tag','shift_t');
      uicontrol('position',[0.03 0.37 0.12 0.05],'tag','compute','enable','off','string','Find ELMs','callback',strcb);
      uicontrol('position',[0.03 0.25 0.12 0.05],'style','edit','tag','filename','enable','off');
      uicontrol('position',[0.03 0.30 0.12 0.05],'tag','savefile','string','Save ASCII','enable','off','callback',strcb);
      uicontrol('position',[0.03 0.18 0.12 0.05],'callback',strcb,'tag','add','string','Add point',...
                 'enable','off','callback',strcb);
      uicontrol('position',[0.03 0.11 0.12 0.05],'callback',strcb,'tag','erase','string','Remove point',...
                 'enable','off','callback',strcb);
      uicontrol('position',[0.03 0.04 0.12 0.05],'tag','multipt','string','   Multi points','style','checkbox','enable','off');
      uicontrol('position',[0.60 0.07 0.35 0.03],'style','text','tag','message','backgroundcolor','g'); 
      uicontrol('position',[0.80 0.92 0.10 0.05],'string','Reset axes','tag','reset','callback',strcb,'enable','off');
      uicontrol('position',[0.650 0.92 0.10 0.05],'string','autoset','tag','autonum','callback',strcb,'value',1,'style','checkbox');
      
      
   case 'download'
   
      s.shot=str2num(get(findobj(gcbf,'tag','shotnumber'),'string'));
      fscope=get(findobj(gcbf,'tag','filterscope'),'string');
      msg=findobj(gcbf,'tag','message');

      if isempty(s.shot)
         set(msg,'string','Please, insert shot number','backgroundcolor','r');
         return;
      end
      shot=double(round(s.shot));
      if or(s.shot<0,length(s.shot)>1)
         set(msg,'string','Shot number must be a positive integer','backgroundcolor','r');
         return;
      end
      set(msg,'string','Connecting to ATLAS...','backgroundcolor','r');
      pause(eps)
      [shoto,dump] = mdsopen('atlas.gat.com::spectroscopy',s.shot);
      if isempty(shoto)
         set(msg,'string','Invalid shot number','backgroundcolor','r');
	 mdsclose;
	 mdsdisconnect;
	 return
      end 
      set(msg,'string',strcat(['Downloading traces from filterscope ',fscope,'...']),'backgroundcolor','r');
      pause(eps)
      s.u = mdsvalue(strcat(['\',fscope]));
      if isstr(s.u)
         set(msg,'string','Invalid filterscope','backgroundcolor','r');
	 mdsclose;
	 mdsdisconnect;
	 return
      end
      s.t = mdsvalue(strcat(['dim_of(\',fscope,')']));
      s.t = s.t/1000;
      s.u = s.u/1e15;
      mdsclose;
      mdsdisconnect;
      
      set(msg,'string',strcat(['Downloaded traces from filterscope ',fscope]),'backgroundcolor','g');
      set(findobj(gcbf,'tag','compute'),'userdata',s,'enable','on');
      set(findobj(gcbf,'tag','sfac'),'enable','on');
      set(findobj(gcbf,'tag','sfac_t'),'enable','on');
      set(findobj(gcbf,'tag','tmax_t'),'enable','on');
      set(findobj(gcbf,'tag','tmax'),'enable','on');
      set(findobj(gcbf,'tag','tmin'),'enable','on');
      set(findobj(gcbf,'tag','tmin_t'),'enable','on');
      set(findobj(gcbf,'tag','shift'),'enable','on');
      set(findobj(gcbf,'tag','shift_t'),'enable','on');
      set(findobj(gcbf,'tag','add'),'enable','off');
      set(findobj(gcbf,'tag','erase'),'enable','off');
      set(findobj(gcbf,'tag','multipt'),'enable','off');
      set(findobj(gcbf,'tag','savefile'),'enable','off');
      set(findobj(gcbf,'tag','filename'),'enable','off');
      set(findobj(gcbf,'tag','reset'),'enable','off')
      delete(findobj(gcbf,'type','axes'));
      
   case 'compute'
   
   %The theoretical factor between the actual value of E(max(du)) and the theoretical formula sigma*sqrt(2*log(N)) is
   %decently approximated by E(max(du))/sigma*sqrt(2*log(N))~=0.92*(1-exp(-(N/0.05).^(0.107))) for 10<N<1e6
   
      sfac=str2num(get(findobj(gcbf,'tag','sfac'),'string'));
      tmin=str2num(get(findobj(gcbf,'tag','tmin'),'string'));
      tmax=str2num(get(findobj(gcbf,'tag','tmax'),'string'));
      shift=str2num(get(findobj(gcbf,'tag','shift'),'string'));
      msg=findobj(gcbf,'tag','message');
      variab={'sfac','tmin','tmax','shift'};
      default={'0.8','0','4','0.1'};
      for i=1:length(variab)
         if isempty(eval(variab{i}))
	    eval(strcat([variab{i},'=',default{i},';']));
	    eval(strcat(['set(findobj(gcbf,''tag'',''',variab{i},'''),''string'',''',default{i},''')']))
	 end
      end	
      an=get(findobj(gcbf,'tag','autonum'),'value');
      
      %converting from milli seconds to seconds 
      shift=shift/1000;
      if tmax<=tmin
         set(msg,'string','t_max must be greater than t_min','backgroundcolor','r');
	 return
      end
      if sfac<0
         set(msg,'string','Statistical factor s-fac should be positive','backgroundcolor','r');
	 return
      end
      if shift<0
         set(msg,'string','Time shift should be positive','backgroundcolor','r');
	 return
      end
      
      s=get(gcbo,'userdata');
      s.ind = find(and(s.t>=tmin,s.t<=tmax));
      u = s.u(s.ind);
      t = s.t(s.ind);
      
      %Implementing first derivative to get rid of secular trend
      du = 0.5*(u(3:end)-u(1:end-2));
      t = t(2:end-1);
      u = u(2:end-1);

      %Smoothing first derivative over 5 points to minimize spurious detections
      %du = conv(du,hann(5));
      %du = du(3:end-2);

      %setting parameters
      vdu = std(du);
      M = length(du);
      
      %Initial step
      facdu = vdu*sqrt(2*log(M))*(an*0.92*(1-exp(-(M/0.05).^(0.107)))+(1-an)*sfac);
      facdu_old = 0;
      ind1 = [];
      i = 0;
      %Iterating to get rid of points above threshold and converge on final variance
      while and(i<100,facdu~=facdu_old)
         
         i = i+1;
         facdu_old = facdu;
         dump2 = find(abs(du)>facdu);
	 dump1 = find(du>facdu);%avoids spurious detections in negative time derivatives
         ind1 = [ind1;dump1];
         du(dump2) = nan;
	 dump=find(~isnan(du));
         vdu = std(du(dump));
	 M=length(dump);
	 %The last term allows sfac to change as M decreases or to keep it fixed
         facdu = vdu*sqrt(2*log(M))*(an*0.92*(1-exp(-(M/0.05).^(0.107)))+(1-an)*sfac);
	 
     end
               
     if isempty(ind1)
        set(msg,'string','No ELM detected, try lowering s-fac','backgroundcolor','r')
	return
     end
     %sorting indices
     ind1 = sort(ind1);
     %Selecting only the first index, i.e. when the ELM is triggered
     %Additionally, selecting ELM at least 0.2 ms apart...to avoid spurious detections
     shift=ceil(shift/mean(diff(t)));
     ind1 = [ind1(1);ind1(find(diff(ind1)>shift)+1)];
     hold off
     set(msg,'string',strcat(['ELMS detected after ',num2str(i),' iterations']),'backgroundcolor','g')
     plot(t,u,'k');
     hold on
     set(gca,'position',[0.25 0.2 0.7 0.7])
     plot(t(ind1),u(ind1),'marker','*','color','r','markersize',14,'linestyle','none')
     legend('Signal','ELM')
     set(gca,'xlim',[tmin-0.001 tmax+0.001])
               
     s.ind1 = ind1;
          
     set(findobj(gcbf,'tag','filename'),'enable','on');
     set(findobj(gcbf,'tag','savefile'),'enable','on','userdata',s); 
     set(findobj(gcbf,'tag','add'),'enable','on','userdata',s);
     set(findobj(gcbf,'tag','erase'),'enable','on','userdata',s); 
     set(findobj(gcbf,'tag','multipt'),'enable','on');
     set(findobj(gcbf,'tag','reset'),'enable','on');
     set(findobj(gcbf,'tag','sfac'),'string',num2str(an*0.92*(1-exp(-(M/0.05).^(0.107)))+(1-an)*sfac));
     
  case 'savefile'
     
     format long
     filename=get(findobj(gcbf,'tag','filename'),'string');
     s=get(gcbo,'userdata');
     if isempty(filename)
        tmin=get(findobj(gcbf,'tag','tmin'),'string');
        tmax=get(findobj(gcbf,'tag','tmax'),'string');
	filename=strcat(['ELMtimes_',num2str(s.shot),'_',tmin,'_',tmax,'.txt']);
     end
     fid=fopen(filename,'w');
     fprintf(fid,'%2.10g\n',s.t(s.ind(s.ind1(:))));
     fclose(fid);
     set(findobj(gcbf,'tag','message'),'string',strcat(['Data saved as ',filename]),'Backgroundcolor','g')
     
  case 'add'
  
     s=get(gcbo,'userdata');
     multipt=get(findobj(gcbf,'tag','multipt'),'value');
     button=1;
     while floor(1/button)
        [a,dump,button]=ginput(1);
	if isempty(a)
	   break
	end
	button=button-multipt+1;%trick to keep the right logic, due to lack of do-while loop in matlab
        [dump,bb]=min(abs(s.t(s.ind)-a));
        s.ind1=[s.ind1;bb];
        s.ind1=sort(s.ind1);
        set(gcbo,'userdata',s);
        set(findobj(gcbf,'tag','erase'),'userdata',s);
        set(findobj(gcbf,'tag','compute'),'userdata',s);
        set(findobj(gcbf,'tag','savefile'),'userdata',s);
        AX=findobj(gcbf,'type','axes');
        AX=get(AX(2),'children');
        set(AX(1),'xdata',s.t(s.ind(s.ind1)),'ydata',s.u(s.ind(s.ind1)));
     end
  
  case 'erase'
  
     s=get(gcbo,'userdata');
     multipt=get(findobj(gcbf,'tag','multipt'),'value');
     button=1;
     while floor(1/button)
        [a,dump,button]=ginput(1);
	if isempty(a)
	   break
	end
	button=button-multipt+1;%trick to keep the right logic, due to lack of do-while loop in matlab
        [dump,bb]=min(abs(s.t(s.ind(s.ind1))-a));
        s.ind1=s.ind1(find([1:length(s.ind1)]-bb));
        set(gcbo,'userdata',s);
        set(findobj(gcbf,'tag','add'),'userdata',s);
        set(findobj(gcbf,'tag','compute'),'userdata',s);
        set(findobj(gcbf,'tag','savefile'),'userdata',s);
        AX=findobj(gcbf,'type','axes');
        AX=get(AX(2),'children');
        set(AX(1),'xdata',s.t(s.ind(s.ind1)),'ydata',s.u(s.ind(s.ind1)));
     end
     
  case 'reset'
  
     AX=findobj(gcf,'type','axes');	
     AX=get(AX(2),'children');
     set(gca,'xlim',[min(get(AX(2),'xdata'))-0.001 max(get(AX(2),'xdata'))+0.001],'ylim',[0 max(get(AX(2),'ydata'))*1.05])
       
end
