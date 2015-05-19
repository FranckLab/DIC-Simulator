function varargout = DIC_evaluator(varargin)
%This program is intended to evaluate the quality of a given DIC image by
%using simulated deformation fields and analyzing the computed
%displacements given the known deformations and speckles of the
%particular-user loaded image
%This program was developed by C. Franck (Brown University), 02/2012 
% DIC_EVALUATOR M-file for DIC_evaluator.fig
%      DIC_EVALUATOR, by itself, creates a new DIC_EVALUATOR or raises the existing
%      singleton*.
%
%      H = DIC_EVALUATOR returns the handle to a new DIC_EVALUATOR or the handle to
%      the existing singleton*.
%
%      DIC_EVALUATOR('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DIC_EVALUATOR.M with the given input arguments.
%
%      DIC_EVALUATOR('Property','Value',...) creates a new DIC_EVALUATOR or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before DIC_evaluator_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to DIC_evaluator_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help DIC_evaluator

% Last Modified by GUIDE v2.5 27-Feb-2012 23:46:29

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @DIC_evaluator_OpeningFcn, ...
                   'gui_OutputFcn',  @DIC_evaluator_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before DIC_evaluator is made visible.
function DIC_evaluator_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to DIC_evaluator (see VARARGIN)

% Choose default command line output for DIC_evaluator
handles.output = hObject;

% define everything here
%handles.input = varargin{1};   %this is the first input into the function
%handles.input2 = varargin{2};
% [x,y] = meshgrid(1:32:512,1:32:512);
% handles.x = x;
% handles.y = y;
% handles.displ = 0.1*x - 0.05*y;
% handles.u = 0.1*x;
% handles.v = -0.05*y;
% 
% 
handles.dicw = 64;
handles.dicd = 16;
handles.trans = 7;
handles.epp = 0.02;


% Update handles structure
guidata(hObject, handles);



% UIWAIT makes DIC_evaluator wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = DIC_evaluator_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in loadimage.
function loadimage_Callback(hObject, eventdata, handles)
% hObject    handle to loadimage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename, pathname] = uigetfile({'*.tif';'*.png';'*.mat';'*.jpg'},'Select your DIC Image')
% if isequal(filename,0)
%    disp('User selected Cancel')
% else
%    disp(['User selected', fullfile(pathname, filename)])
% end
%  dvc_file = load(filename);
%  handles.pic = dvc_file.vol_stack; 
handles.pic = importdata([pathname,filename]);
axes(handles.axes1);
imagesc(handles.pic);
set(gca,'fontweight','b','FontSize',12,'xcolor','w','ycolor','w');
[m,n,p] = size(handles.pic);
if p ~= 1,
    handles.pic = squeeze(sum(handles.pic(:,:,1:p),3));
end
guidata(hObject,handles);


function dicw_Callback(hObject, eventdata, handles)
% hObject    handle to dicw (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dicw as text
%        str2double(get(hObject,'String')) returns contents of dicw as a double
w = str2double(get(hObject, 'String'));
if isnan(w)
    set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
end

% Save the new density value
handles.dicw = w;
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function dicw_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dicw (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
handles.dicw = 64;
guidata(hObject,handles);



function dicd_Callback(hObject, eventdata, handles)
% hObject    handle to dicd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dicd as text
%        str2double(get(hObject,'String')) returns contents of dicd as a double
d = str2double(get(hObject, 'String'));
if isnan(d)
    set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
end

% Save the new density value
handles.dicd = d;
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function dicd_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dicd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
handles.dicd = 16;
guidata(hObject,handles);



function epp_Callback(hObject, eventdata, handles)
% hObject    handle to epp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of epp as text
%        str2double(get(hObject,'String')) returns contents of epp as a double
epp = str2double(get(hObject, 'String'));
if isnan(epp)
    set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
end

% Save the new density value
handles.epp = epp;
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function epp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to epp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
handles.epp = 0.02;
guidata(hObject,handles);



function trans_Callback(hObject, eventdata, handles)
% hObject    handle to trans (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of trans as text
%        str2double(get(hObject,'String')) returns contents of trans as a double
trans = str2double(get(hObject, 'String'));
if isnan(trans)
    set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
end

% Save the new density value
handles.trans = trans;
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function trans_CreateFcn(hObject, eventdata, handles)
% hObject    handle to trans (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
handles.trans = 7;
guidata(hObject,handles);

% --- Executes on button press in gobutton.
function gobutton_Callback(hObject, eventdata, handles)
% hObject    handle to gobutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% axes(handles.axes1);
% handles.x;
% handles.y;
% handles.displ;
% contourf(handles.x, handles.y,handles.displ);

%Simulate deformed speckle images
handles.epp;
[handles.tpic handles.strpic handles.shpic handles.ppic dispp] = speckles(handles.epp, handles.trans, handles.pic);

%Execute 2D DIC code *******************************************************
%handles.inc=0;
%Translation
[ut vt] = run2dic(handles.dicw,handles.dicd,handles.pic,handles.tpic);
[dudxt, dudyt, dvdxt, dvdyt]=grad2d(ut,vt,1,1,handles.dicd*0.35,handles.dicd*0.35);
handles.ut = ut; 
handles.vt = vt;
handles.e11t = dudxt; 
handles.e22t = dvdyt;
handles.e12t = dudyt; 
handles.e21t = dvdxt;

%Uniaxal Tension
[ua va] = run2dic(handles.dicw,handles.dicd,handles.pic,handles.strpic);
[dudxa, dudya, dvdxa, dvdya]=grad2d(ua,va,1,1,handles.dicd*0.35,handles.dicd*0.35);
handles.ua = ua; 
handles.va = va;
handles.e11a = dudxa; 
handles.e22a = dvdya;
handles.e12a = dudya; 
handles.e21a = dvdxa;

%Pure Shear
[us vs] = run2dic(handles.dicw,handles.dicd,handles.pic,handles.shpic);
[dudxs, dudys, dvdxs, dvdys]=grad2d(us,vs,1,1,handles.dicd*0.35,handles.dicd*0.35);
handles.us = us; 
handles.vs = vs;
handles.e11s = dudxs; 
handles.e22s = dvdys;
handles.e12s = dudys; 
handles.e21s = dvdxs;

%Point Force
[up vp] = run2dic(handles.dicw,handles.dicd,handles.pic,handles.ppic);
[dudxp, dudyp, dvdxp, dvdyp]=grad2d(up,vp,1,1,handles.dicd*0.35,handles.dicd*0.35);
handles.up = up; 
handles.vp = vp;
handles.e11p = dudxp; 
handles.e22p = dvdyp;
handles.e12p = dudyp; 
handles.e21p = dvdxp;

[M,N]=size(handles.pic);%Find the size of the confocal image
%w0=64;d0=16; %w0:Subset size & d0:Subset spacing, 
    %For example, when M=N=P=512, (w0=64,d0=32) gives 15x15x15 data points for correlations
x0=N/2;y0=M/2; %find the center of the image
d0=handles.dicd;
%y0=256;x0=256;z0=256; %Manually select center
m=round(M/2/d0)-1;n=round(N/2/d0)-1;% number of correlation points in half of the image (any integer)
[x,y]=meshgrid([-n:n]*d0+x0,[-m:m]*d0+y0);%
handles.x = x;
handles.y = y;

if handles.popup == 1,
    handles.u = handles.ut;
    handles.v = handles.vt;
    handles.e11 = handles.e11t; 
    handles.e22 = handles.e22t;
    handles.e12 = handles.e12t; 
    handles.e21 = handles.e21t;
end
if handles.popup == 2,
    handles.u = handles.ua;
    handles.v = handles.va;
    handles.e11 = handles.e11a; 
    handles.e22 = handles.e22a;
    handles.e12 = handles.e12a; 
    handles.e21 = handles.e21a;
end
if handles.popup == 3,
    handles.u = handles.us;
    handles.v = handles.vs;
    handles.e11 = handles.e11s; 
    handles.e22 = handles.e22s;
    handles.e12 = handles.e12s; 
    handles.e21 = handles.e21s;
end

handles.displ = sqrt(handles.u.^2 + handles.v.^2);

axes(handles.axes2);
%im(handles.tpic);
handles.x;
handles.y;
%handles.displ;
contourf(handles.x, handles.y,handles.u); colorbar;
h = colorbar;
set(get(h,'title'),'string','Pixel','fontweight','b','fontsize',12,'color','w');
set(h,'fontweight','b', 'fontsize',12);
grid on
axis tight
set(gca,'fontweight','b','FontSize',12,'xcolor','w','ycolor','w');


axes(handles.axes3);
handles.x;
handles.y;
%handles.displ;
contourf(handles.x, handles.y,handles.v); colorbar;
h = colorbar;
set(get(h,'title'),'string','Pixel','fontweight','b','fontsize',12,'color','w');
set(h,'fontweight','b', 'fontsize',12);
grid on
axis tight
set(gca,'fontweight','b','FontSize',12,'xcolor','w','ycolor','w');


axes(handles.axes5);
clear err;
err = (handles.ut - dispp.ut2)./dispp.ut2*100;
err = err(:);
hist(err,20);
set(gca,'fontweight','b','FontSize',10,'xcolor','k','ycolor','k');
%xlabel('difference (pixel)','fontweight','b','FontSize',10);
ylabel('U_1','fontweight','b','FontSize',10);

axes(handles.axes6);
%Construct the solution right here:
[M,N]=size(handles.tpic);
x0=N/2;y0=M/2; %find the center of the image
d0 = handles.dicd;
m=round(M/2/d0)-1;n=round(N/2/d0)-1;% number of correlation points in half of the image (any integer)
[x,y]=meshgrid([-n:n]*d0+x0,[-m:m]*d0+y0);
strain = handles.epp;
ua = strain*x;   % apply strain along x-direction
va = -strain/2*y;   % Poisson's contraction with nu = 0.5
err = (handles.ua - ua)./ua*100;
err = err(:);
hist(err,20);
set(gca,'fontweight','b','FontSize',10,'xcolor','k','ycolor','k');
%xlabel('difference (pixel)','fontweight','b','FontSize',10);
%ylabel('Frequency','fontweight','b','FontSize',10);

axes(handles.axes7);
us = strain*y;   %  shear in x
vs = strain*x;   %  shear in y
err = (handles.us - us)./us*100;
err = err(:);
hist(err,20);
set(gca,'fontweight','b','FontSize',10,'xcolor','k','ycolor','k');
%xlabel('difference (pixel)','fontweight','b','FontSize',10);
%ylabel('Frequency','fontweight','b','FontSize',10);

axes(handles.axes8);
% generate point force image
nu = 0.45;
E = 1000;
Fx = E*strain;
[m,n] = size(x);

x = x-round(N/2);
y = y - round(M/2);
r = sqrt(x.^2 + y.^2);
ind = find(r == 0); %Exclude the singularity at r = 0
r(ind) = 5*10^-4;

coeff = (1+nu)/(2*pi*E);
z = 0; Fz = 0; Fy = Fx;
up = coeff*( (x.*z./r.^3 - (1-2*nu)*x./( r.*(r + z)))*Fz...
    + (2*(1-nu)*r + z)./(r.*(r + z))*Fx...
    + (2*r.*(nu*r + z) + z.^2).*x./(r.^3.*(r + z).^2).*(x*Fx + y*Fy));

vp = coeff*( (y.*z./r.^3 - (1-2*nu)*y./( r.*(r + z)))*Fz...
    + (2*(1-nu)*r + z)./(r.*(r + z))*Fy...
    + (2*r.*(nu*r + z) + z.^2).*y./(r.^3.*(r + z).^2).*(x*Fx + y*Fy));
err = (handles.up - up)./up*100;
err = err(:);
hist(err,20);
set(gca,'fontweight','b','FontSize',10,'xcolor','k','ycolor','k');
%xlabel('difference (pixel)','fontweight','b','FontSize',10);
%ylabel('Frequency','fontweight','b','FontSize',10);

axes(handles.axes9);
if dispp.vt2 == 0,
    err = (handles.vt);
else
err = (handles.vt - dispp.vt2)./dispp.vt2*100;
end
err = err(:);
hist(err,20);
set(gca,'fontweight','b','FontSize',10,'xcolor','k','ycolor','k');
%xlabel('difference (pixel)','fontweight','b','FontSize',10);
ylabel('U_2','fontweight','b','FontSize',10);

axes(handles.axes10);
err = (handles.va - va)./va*100;
err = err(:);
hist(err,20);
set(gca,'fontweight','b','FontSize',10,'xcolor','k','ycolor','k');
%xlabel('difference (pixel)','fontweight','b','FontSize',10);
%ylabel('Frequency','fontweight','b','FontSize',10);

axes(handles.axes11);
err = (handles.vs - vs)./vs*100;
err = err(:);
hist(err,20);
set(gca,'fontweight','b','FontSize',10,'xcolor','k','ycolor','k');
%xlabel('difference (pixel)','fontweight','b','FontSize',10);
%ylabel('Frequency','fontweight','b','FontSize',10);

axes(handles.axes12);
err = (handles.vp - vp)./vp*100;
err = err(:);
hist(err,20);
set(gca,'fontweight','b','FontSize',10,'xcolor','k','ycolor','k');
%xlabel('difference (pixel)','fontweight','b','FontSize',10);
%ylabel('Frequency','fontweight','b','FontSize',10);


% handles.x;
% handles.x;
% handles.x;
% handles.x;
% handles.x;
% handles.y;
% handles.u;
% handles.v;
% quiver(handles.x, handles.y,handles.u,handles.v);
% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in magu.
function magu_Callback(hObject, eventdata, handles)
% hObject    handle to magu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.popup == 1,
    handles.u = handles.ut;
    handles.v = handles.vt;
end
if handles.popup == 2,
    handles.u = handles.ua;
    handles.v = handles.va;
end
if handles.popup == 3,
    handles.u = handles.us;
    handles.v = handles.vs;
end
if handles.popup == 4,
    handles.u = handles.up;
    handles.v = handles.vp;
end
handles.displ = sqrt(handles.u.^2 + handles.v.^2);

axes(handles.axes2);
cla;
contourf(handles.x, handles.y,handles.displ); 
h = colorbar;
set(get(h,'title'),'string','Pixel','fontweight','b','fontsize',12,'color','w');
set(h,'fontweight','b', 'fontsize',12);
grid on
axis tight
set(gca,'fontweight','b','FontSize',12,'xcolor','w','ycolor','w');

axes(handles.axes3);
cla;
contourf(handles.x, handles.y,handles.displ); 
hold on;
quiver(handles.x, handles.y,handles.u,handles.v,'w'); 
h = colorbar;
set(get(h,'title'),'string','Pixel','fontweight','b','fontsize',12,'color','w');
set(h,'fontweight','b', 'fontsize',12);
grid on
axis tight
set(gca,'fontweight','b','FontSize',12,'xcolor','w','ycolor','w');
% Hint: get(hObject,'Value') returns toggle state of magu


% --- Executes on button press in ucomp.
function ucomp_Callback(hObject, eventdata, handles)
% hObject    handle to ucomp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.popup == 1,
    handles.u = handles.ut;
    handles.v = handles.vt;
end
if handles.popup == 2,
    handles.u = handles.ua;
    handles.v = handles.va;
end
if handles.popup == 3,
    handles.u = handles.us;
    handles.v = handles.vs;
end
if handles.popup == 4,
    handles.u = handles.up;
    handles.v = handles.vp;
end

axes(handles.axes2);
cla;
contourf(handles.x, handles.y,handles.u); 
h = colorbar;
set(get(h,'title'),'string','Pixel','fontweight','b','fontsize',12,'color','w');
set(h,'fontweight','b', 'fontsize',12);
grid on
axis tight
set(gca,'fontweight','b','FontSize',12,'xcolor','w','ycolor','w');

axes(handles.axes3);
cla;
contourf(handles.x, handles.y,handles.v); 
h = colorbar;
set(get(h,'title'),'string','Pixel','fontweight','b','fontsize',12,'color','w');
set(h,'fontweight','b', 'fontsize',12);
grid on
axis tight
set(gca,'fontweight','b','FontSize',12,'xcolor','w','ycolor','w');
% Hint: get(hObject,'Value') returns toggle state of ucomp

%------------------------------------------------------
function [tpic strpic shpic ppic dispp] = speckles(strain, trans, pic);

 %General image warping based on deformation field                   
C= pic; % test image
[x, y] = meshgrid(1:size(C,2), 1:size(C,1));

% generate translated image
u = trans;   % x-translation
v = 0;   % y-translation
dispp.ut2 = u;
dispp.vt2 = v;
tpic = interp2(double(C), x-u, y-v);
%figure,im(tpic);
clear u; clear v;

% generate uniaxially stretched image
u = strain*x;   % apply strain along x-direction
v = -strain/2*y;   % Poisson's contraction with nu = 0.5
dispp.ua2 = u;
dispp.va2 = v;
strpic = interp2(double(C), x-u, y-v);
%figure,im(strpic);
clear u; clear v;

% generate pure sheared  image
u = strain*y;   %  shear in x
v = strain*x;   %  shear in y
dispp.us2 = u;
dispp.vs2 = v;
% compute the warped image - the subtractions are because we're specifying
% where in the original image each pixel in the new image comes from
shpic = interp2(double(C), x-u, y-v);
%figure,im(shpic);
clear u; clear v;

% generate point force image
nu = 0.45;
E = 1000;
Fx = E*strain;
[m,n] = size(x);

x = x-round(size(C,2)/2);
y = y - round(size(C,1)/2);
r = sqrt(x.^2 + y.^2);
ind = find(r == 0); %Exclude the singularity at r = 0
r(ind) = 5*10^-4;

coeff = (1+nu)/(2*pi*E);
z = 0; Fz = 0; Fy = Fx;
u = coeff*( (x.*z./r.^3 - (1-2*nu)*x./( r.*(r + z)))*Fz...
    + (2*(1-nu)*r + z)./(r.*(r + z))*Fx...
    + (2*r.*(nu*r + z) + z.^2).*x./(r.^3.*(r + z).^2).*(x*Fx + y*Fy));

v = coeff*( (y.*z./r.^3 - (1-2*nu)*y./( r.*(r + z)))*Fz...
    + (2*(1-nu)*r + z)./(r.*(r + z))*Fy...
    + (2*r.*(nu*r + z) + z.^2).*y./(r.^3.*(r + z).^2).*(x*Fx + y*Fy));
max(u(:));
dispp.up2 = u;
dispp.vp2 = v;
% compute the warped image - the subtractions are because we're specifying
% where in the original image each pixel in the new image comes from
ppic = interp2(x,y,double(C), x-u, y-v);
%figure,im(pic);
%figure,im(ppic);
v;



%----------------------------------------------------------

% --- Executes on selection change in output_selection.
function output_selection_Callback(hObject, eventdata, handles)
% hObject    handle to output_selection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns output_selection contents as cell array
%        contents{get(hObject,'Value')} returns selected item from output_selection
popup_sel_index = get(handles.output_selection, 'Value');
handles.popup = popup_sel_index;
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function output_selection_CreateFcn(hObject, eventdata, handles)
% hObject    handle to output_selection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
handles.popup = 1;
guidata(hObject, handles);

% --- Executes on button press in diagstrains.
function diagstrains_Callback(hObject, eventdata, handles)
% hObject    handle to diagstrains (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.popup == 1,
    handles.e11 = handles.e11t; 
    handles.e22 = handles.e22t;
end
if handles.popup == 2,
    handles.e11 = handles.e11a; 
    handles.e22 = handles.e22a;
end
if handles.popup == 3,
    handles.e11 = handles.e11s; 
    handles.e22 = handles.e22s;
end
if handles.popup == 4,
    handles.e11 = handles.e11p; 
    handles.e22 = handles.e22p;
end
% Hint: get(hObject,'Value') returns toggle state of diagstrains
axes(handles.axes2);
cla;
contourf(handles.x, handles.y,handles.e11); 
h = colorbar;
set(get(h,'title'),'string','Pixel','fontweight','b','fontsize',12,'color','w');
set(h,'fontweight','b', 'fontsize',12);
grid on
axis tight
set(gca,'fontweight','b','FontSize',12,'xcolor','w','ycolor','w');

axes(handles.axes3);
cla;
contourf(handles.x, handles.y,handles.e22); 
h = colorbar;
set(get(h,'title'),'string','Pixel','fontweight','b','fontsize',12,'color','w');
set(h,'fontweight','b', 'fontsize',12);
grid on
axis tight
set(gca,'fontweight','b','FontSize',12,'xcolor','w','ycolor','w');

% --- Executes on button press in divstrains.
function divstrains_Callback(hObject, eventdata, handles)
% hObject    handle to divstrains (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of divstrains
if handles.popup == 1,
    handles.e12 = handles.e12t; 
    handles.e21 = handles.e21t;
end
if handles.popup == 2,
    handles.e12 = handles.e12a; 
    handles.e21 = handles.e21a;
end
if handles.popup == 3,
    handles.e12 = handles.e12s; 
    handles.e21 = handles.e21s;
end
if handles.popup == 4,
    handles.e12 = handles.e12p; 
    handles.e21 = handles.e21p;
end

axes(handles.axes2);
cla;
contourf(handles.x, handles.y,handles.e12); 
h = colorbar;
set(get(h,'title'),'string','Pixel','fontweight','b','fontsize',12,'color','w');
set(h,'fontweight','b', 'fontsize',12);
grid on
axis tight
set(gca,'fontweight','b','FontSize',12,'xcolor','w','ycolor','w');

axes(handles.axes3);
cla;
contourf(handles.x, handles.y,handles.e21); 
h = colorbar;
set(get(h,'title'),'string','Pixel','fontweight','b','fontsize',12,'color','w');
set(h,'fontweight','b', 'fontsize',12);
grid on
axis tight
set(gca,'fontweight','b','FontSize',12,'xcolor','w','ycolor','w');

%*****************************************************************RUN DIC

%This is a 2D correlation code derived from the 3D code

function [ud vd] = run2dic(w0,d0,ia,ib);

inc = 0; stack = 1;

%--------------------------Construct 2d grid points for measurements-----------
[M,N]=size(ia);%Find the size of the confocal image
%w0=64;d0=16; %w0:Subset size & d0:Subset spacing, 
    %For example, when M=N=P=512, (w0=64,d0=32) gives 15x15x15 data points for correlations
x0=N/2;y0=M/2; %find the center of the image
%y0=256;x0=256;z0=256; %Manually select center
m=round(M/2/d0)-1;n=round(N/2/d0)-1;% number of correlation points in half of the image (any integer)
[x,y]=meshgrid([-n:n]*d0+x0,[-m:m]*d0+y0);%construct a set of 3d grid points centered at the center of the confocal image. 
ua=zeros(size(x));va=ua;%initialize displacement variables, size(x) is odd number
 [ud,vd]=dvc_full_2d(ia,ib,x,y,ua,va,w0,inc);
 
 %--------FFT-based correlation ---------------------------
 
 function [u,v]=dvc_full_2d(m0,m1,x0,y0,ua,va,w0,inc)
[M N]=size(m0);
[m n]=size(x0);
u=zeros(m,n)*nan;v=u;
%u2=zeros(m,n)*nan;v2=u2;
x0=round(x0);y0=round(y0);
x1=x0+round(ua);y1=y0+round(va);
if inc == 1,
    x0 = x1; y0 = y1; 
    ua=zeros(size(x0));va=ua;
end
h = waitbar(0,'Running 2D DIC Code... ');
for i=1:m,
   waitbar(i/m,h);
   for j=1:n,
      i0=(y0(i,j)-w0/2+1); i1=y0(i,j)+w0/2;
      j0=(x0(i,j)-w0/2+1); j1=x0(i,j)+w0/2;
      p0=(y1(i,j)-w0/2+1); p1=y1(i,j)+w0/2;
      q0=(x1(i,j)-w0/2+1); q1=x1(i,j)+w0/2;
      if i0>=1 & i0<=M & j0>=1 & j0<=N & ...
             i1>=1 & i1<=M & j1>=1 & j1<=N & ...
             p0>=1 & p0<=M & q0>=1 & q0<=N & ...
             p1>=1 & p1<=M & q1>=1 & q1<=N,
             dm0=m0(i0:i1,j0:j1);
             dm1=m1(p0:p1,q0:q1);
             [du dv dc]=decorr2a(dm0,dm1);
             %[du,dv] = normxcorr2D(dm0,dm1);
             if(length(du)*length(dv)==1)
%                x(i,j,k)=(j-1)*d+d/2;
%                y(i,j,k)=(i-1)*d+d/2;
%                z(i,j,k)=(k-1)*d+d/2;
                u(i,j)=round(ua(i,j))+du;
                v(i,j)=round(va(i,j))+dv;
             end
       end
   end
end
close(h);

%-----------------------------------
function [dudx, dudy, dvdx, dvdy]=grad2d(u,v,ucal,vcal,xcal,ycal)
[M N]=size(u);

[dx,dy]=meshgrid(-[-1:1]/18,-[-1:1]/18);
dudx=ones(M,N)*nan;
dudy=dudx;
dvdx=dudx;
dvdy=dudx;

dudx(2:M-1,2:N-1)=convn(u,dx,'valid')*ucal/xcal;
dudy(2:M-1,2:N-1)=convn(u,dy,'valid')*ucal/ycal;

dvdx(2:M-1,2:N-1)=convn(v,dx,'valid')*vcal/xcal;
dvdy(2:M-1,2:N-1)=convn(v,dy,'valid')*vcal/ycal;

%-----------------------------------

function [du,dv,dc]=decorr2a(dm0,dm1)
[M,N]=size(dm0);
dm0=double(dm0); dm1=double(dm1);
c=abs(fftshift(ifft2(conj(fft2(dm0)).*fft2(dm1))));
cn=abs(fftshift(ifft2(conj(fft2(dm0))./abs(fft2(dm0)).*fft2(dm1)./abs(fft2(dm1)))));
   cc=sum(sum((dm0-dm1).*(dm0-dm1)));
   if cc==0 | isnan(sum(c(:))),
   	du=nan;
   	dv=nan;
    dc=nan;
   else 
      [dv,du]=find(c==max(c(:)));

      if dv==1|dv==M|du==1|du==N,
            	du=nan;
                dv=nan;
                dc=nan;
      else
        cc=c(dv-1:dv+1,du-1:du+1);
        [mm,nn]=size(cc);
        [xx,yy]=meshgrid(-(mm-1)/2:(mm-1)/2,-(nn-1)/2:(nn-1)/2);
        xx2=xx.*xx;
        yy2=yy.*yy;
        A= [sum(sum(xx2.*xx2)) sum(sum(xx2.*yy2)) 0 0 sum(sum(xx2));
            sum(sum(xx2.*yy2)) sum(sum(yy2.*yy2)) 0 0 sum(sum(yy2));
            0 0 sum(sum(xx2)) 0 0;
            0 0 0 sum(sum(yy2)) 0;
            sum(sum(xx2)) sum(sum(yy2)) 0 0 sum(sum(ones(mm,nn)))];
        y=[sum(sum(xx2.*cc)) sum(sum(yy2.*cc)) sum(sum(xx.*cc)) sum(sum(yy.*cc)) sum(sum(cc))]';
        x=A\y;
        x0=-x(3)/2/x(1);
        y0=-x(4)/2/x(2);
        dc=sqrt(x(1)^2+x(2)^2);

      du=du+x0-(N/2+1);
      dv=dv+y0-(M/2+1);
    end 
  end


 





