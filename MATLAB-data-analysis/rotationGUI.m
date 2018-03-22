function [xshift, yshift, rot, transparency_red, transparency_green, overlay_ref, overlay_bkg] = rotationGUI(I, I_REF, CELLS_X, CELLS_Y, CELLS_REF_X, CELLS_REF_Y, MIN, MIN_REF, day, exp, transparency_red, transparency_green, overlay_ref, overlay_bkg)
    %# read image
    % I = imread(img);
    
    maxshift = 200;
    filtsize = 3;
    
%     transparency_red = 0.3;
%     transparency_green = 0.6;
%     overlay_ref = 0.35;
%     overlay_bkg = 0.45;
    
    
    radius = 30;%8; % HPF
    amt = 1200;%230; % HPF
    

    width = size(I,2);
    height = size(I,1);
    CELLS_X_ROT = zeros(1,length(CELLS_X));
    CELLS_Y_ROT = zeros(1,length(CELLS_Y));
    CELLS_X_XSH = zeros(1,length(CELLS_X));
    CELLS_Y_YSH = zeros(1,length(CELLS_Y));
   
 
    xshift = 0;
    yshift = 0;
    rot = 0;
    
    
    GAUSS = 20*gauss_filt(filtsize,1);
%     MIN = double(imfilter(MIN, GAUSS, 'same'));
%     MIN_REF = double(imfilter(MIN_REF, GAUSS, 'same'));
    
    MIN = double(MIN);
    MIN_REF = double(MIN_REF);
     
    
    MIN = MIN./max(max(MIN));
    MIN_REF = MIN_REF./max(max(MIN_REF));
    
    
    
    % hh = fspecial('sobel',50,45); 
    % MIN_REF = imfilter(MIN_REFB,hh,'same');
    
   

    
    
    % Enhance background
    
    % H = 1 - fspecial('gaussian' ,[5 5],2); % create unsharp mask
    H = padarray(4,[2 2]) - fspecial('gaussian' ,[5 5],0.5); % create unsharp mask
    % MIN_REF = imfilter(MIN_REF,H);  % create a sharpened version of the image using that mask
    MIN_REF = imsharpen(MIN_REF,'Radius',radius,'Amount',amt);
    MIN_REF = MIN_REF./max(max(MIN_REF));
    
    
    % MIN = imfilter(MIN,H);  % create a sharpened version of the image using that mask
    MIN = imsharpen(MIN,'Radius',radius,'Amount',amt);
    MIN = MIN./max(max(MIN));

    MIN = 1-MIN;
    MIN_REF = 1-MIN_REF;
    
    
    
    
    
    
    
%     Gx = [-1 1];
%     Gy = Gx';
%     Ix = conv2(MIN,Gx,'same');
%     Iy = conv2(MIN,Gy,'same');
%     figure; image(Ix+Iy,'CDataMapping','scaled'); colormap('gray'); title('Ix');
    
     
    [~, threshold] = edge(MIN, 'sobel');
    fudgeFactor = 1.4;%1;
    
    EDGES = edge(MIN,'sobel', threshold * fudgeFactor);
    EDGES_REF = edge(MIN_REF,'sobel', threshold * fudgeFactor);
    
    EDGES(1:filtsize,:) = 0;
    EDGES((end-filtsize):end,:) = 0;
    EDGES(:, 1:filtsize) = 0;
    EDGES(:,end-filtsize:end) = 0;
    

    EDGES_REF(1:filtsize,:) = 0;
    EDGES_REF((end-filtsize):end,:) = 0;
    EDGES_REF(:, 1:filtsize) = 0;
    EDGES_REF(:,end-filtsize:end) = 0;
    
  
    


    hFig = figure('menu','none');
    
    title(strcat('Day ', num2str(day), ', exp ', num2str(exp)));
      
      
    hAx = axes('Parent',hFig);
    set(hFig, 'Position',[50,50,3.5*size(I,2),3.5*size(I,1)]);
    
    red = cat(3, ones(size(I,1), size(I,2)), zeros(size(I,1), size(I,2)), zeros(size(I,1), size(I,2)));
    green = cat(3, zeros(size(I,1), size(I,2)), ones(size(I,1), size(I,2)), zeros(size(I,1), size(I,2)));
    blue = cat(3, zeros(size(I,1), size(I,2)), zeros(size(I,1), size(I,2)), ones(size(I,1), size(I,2)));
    
    black = zeros(size(I,1), size(I,2));


    
    uicontrol('Parent',hFig, 'Style','slider', 'Value',0, 'Min',-180,...
        'Max',180, 'SliderStep',[0.25 5]./360, ...
        'Position',[80 125 3.2*size(I,2) 20], 'Callback',@slider_callback_rot) 
    hTxt_rot = uicontrol('Style','text', 'Position',[5 125 80 15], 'String',strcat('Rotation: ','0'));
     
    uicontrol('Parent',hFig, 'Style','slider', 'Value',0, 'Min',-maxshift+1,...
        'Max',maxshift, 'SliderStep',[1 10]./(2*maxshift), ...
        'Position',[80 105 3.2*size(I,2) 20], 'Callback',@slider_callback_x) 
    hTxt_x = uicontrol('Style','text', 'Position',[5 105 80 15], 'String',strcat('X shift: ','0'));
     
    uicontrol('Parent',hFig, 'Style','slider', 'Value',0, 'Min',-maxshift+1,...
        'Max',maxshift, 'SliderStep',[1 10]./(2*maxshift), ...
        'Position',[80 85 3.2*size(I,2) 20], 'Callback',@slider_callback_y) 
    hTxt_y = uicontrol('Style','text', 'Position',[5 85 80 15], 'String',strcat('Y shift: ','0'));
    

    uicontrol('Parent',hFig, 'Style','slider', 'Value',transparency_red, 'Min',0,...
        'Max',1, 'SliderStep',[1 10]./100, ...
        'Position',[80 65 3.2*size(I,2) 20], 'Callback',@slider_callback_transpr) 
    hTxt_transpr = uicontrol('Style','text', 'Position',[5 65 80 15], 'String',strcat('Red: ',num2str(transparency_red,'%.02f')));
    
    uicontrol('Parent',hFig, 'Style','slider', 'Value',transparency_green, 'Min',0,...
        'Max',1, 'SliderStep',[1 10]./100, ...
        'Position',[80 45 3.2*size(I,2) 20], 'Callback',@slider_callback_transpg) 
    hTxt_transpg = uicontrol('Style','text', 'Position',[5 45 80 15], 'String',strcat('Green: ',num2str(transparency_green,'%.02f')));
    
    uicontrol('Parent',hFig, 'Style','slider', 'Value',overlay_bkg, 'Min',0,...
        'Max',1, 'SliderStep',[1 10]./100, ...
        'Position',[80 25 3.2*size(I,2) 20], 'Callback',@slider_callback_transpb) 
    hTxt_transpb = uicontrol('Style','text', 'Position',[5 25 80 15], 'String',strcat('BCKG: ',num2str(overlay_bkg,'%.02f')));
    
    

    

    
   
    
  

    % imshow(AVG, 'Parent',hAx); hold on;
    hold on;

    hh = imshow(red);
    set(hh, 'AlphaData', transparency_red.*I_REF); 
    
    h = imshow(green);
    set(h, 'AlphaData', transparency_green.*I); 
    
    
    min_ref_h = imshow(blue);
    set(min_ref_h, 'AlphaData', overlay_ref.*MIN_REF); 
    
    
    min_h = imshow(black);
    set(min_h, 'AlphaData', overlay_bkg.*MIN); 
    
    
    % ref_edges = imagesc(1-EDGES_REF); colormap(gray); hold on;
    % set(ref_edges, 'AlphaData', 0.03*ones(size(EDGES_REF,1),size(EDGES_REF,2))); 
    
    % hhh = imshow(blue);
    % set(hhh, 'AlphaData', transparency.*EDGES); 
   
    
%     ref_bkg = imagesc(1-MIN_REF); colormap(gray); hold on;
%     set(ref_bkg, 'AlphaData', 0.6*ones(size(ref_bkg,1),size(ref_bkg,2))); 
%     
%     
%     bkg = imagesc(1-MIN); colormap(gray); hold on;
%     set(bkg, 'AlphaData', 0.6*ones(size(bkg,1),size(bkg,2))); 
    
    
    plot(CELLS_REF_X,CELLS_REF_Y,'.r','Markersize',10);
    plot_xy = plot(CELLS_X,CELLS_Y,'.g');
    
    
    
    
    img = I;
    img2 = I;
    img3 = I;
    img_0 = MIN;
    img_1 = MIN;
    img_2 = MIN;


    
    edge1 = EDGES;
    edge2 = EDGES;
    edge2 = EDGES;
    
    
  
    uiwait(hFig)

    
    
    
    















    function slider_callback_rot(hObj, eventdata)
        img = I;
        img2 = I;
        img3 = I;
        img_0 = MIN;
        img_1 = MIN;
        img_2 = MIN;
        edge1 = EDGES;
        edge2 = EDGES;
        edge2 = EDGES;

        reset(plot_xy);
        
        
        rot = get(hObj,'Value');  % round  
        set(hTxt_rot, 'String',strcat('Rotation: ',num2str(rot,'%.02f'))); 
%         xshift = round(get(hObj,'Value'));
%         set(hTxt_x, 'String', strcat('X shift: ',num2str(xshift,'%.02f')));
%         yshift = round(get(hObj,'Value'));
%         set(hTxt_y, 'String', strcat('Y shift: ',num2str(yshift,'%.02f')));
%         
%         transparency_green = get(hObj,'Value');
%         set(hTxt_transpg, 'String',strcat('Green: ',num2str(transparency_green,'%.02f')));
%         transparency_red = get(hObj,'Value');
%         set(hTxt_transpr, 'String',strcat('Red: ',num2str(transparency_red,'%.02f')));
%         overlay_bkg = get(hObj,'Value');
%         set(hTxt_transpb, 'String',strcat('BCK: ',num2str(overlay_bkg,'%.02f')));
%         overlay_ref = get(hObj,'Value');
%         set(hTxt_transpref, 'String',strcat('BCK: ',num2str(overlay_ref,'%.02f')));

        
        % Rotation
    
        img = imrotate(I,rot,'bilinear','crop'); 
        img_0 = imrotate(MIN,rot,'bilinear','crop');
        edge1 = imrotate(EDGES,rot,'bilinear','crop'); 
        edge2 = edge1;
        img2 = img;
        set(h, 'AlphaData', transparency_green.*img); 
        % set(hhh, 'AlphaData', transparency.*edge1); 
        
        set(min_h, 'AlphaData', overlay_bkg.*img_0);    
        rot_matrix = [cos(-rot*pi/180) -sin(-rot*pi/180); sin(-rot*pi/180) cos(-rot*pi/180)];
        for k = 1:length(CELLS_X)
            [coord] = [width/2; height/2] + rot_matrix*([CELLS_X(k); CELLS_Y(k)]-[width/2; height/2]);
            CELLS_X_ROT(k) = coord(1);
            CELLS_Y_ROT(k) = coord(2);
        end
        
        
        
        % X shift
        
        CELLS_X_XSH = CELLS_X_ROT - xshift;
        if (xshift > 0)
            img2 = zeros(size(img,1),size(img,2));
            img2(:,1:end-xshift) = img(:,xshift+1:end);
            img_1 = zeros(size(img,1),size(img,2));
            img_1(:,1:end-xshift) = img_0(:,xshift+1:end);
            edge2 = zeros(size(img,1),size(img,2));
            edge2(:,1:end-xshift) = edge1(:,xshift+1:end);
        elseif (xshift < 0)
            img2 = zeros(size(img,1),size(img,2));
            img2(:,-xshift+1:end) = img(:,1:end+xshift);
            img_1 = zeros(size(img,1),size(img,2));
            img_1(:,-xshift+1:end) = img_0(:,1:end+xshift);
            edge2 = zeros(size(img,1),size(img,2));
            edge2(:,-xshift+1:end) = edge1(:,1:end+xshift);
        else
            img2 = img;
            edge2 = edge1;
            img_1 = img_0;
        end  
        
        
        
        % Y shift
        
        CELLS_Y_YSH = CELLS_Y_ROT - yshift;
        
        if (yshift > 0)
            img3 = zeros(size(img2,1),size(img2,2));
            img3(1:end-yshift,:) = img2(yshift+1:end,:);
            img_2 = zeros(size(img2,1),size(img2,2));
            img_2(1:end-yshift,:) = img_1(yshift+1:end,:);
            edge3 = zeros(size(img2,1),size(img2,2));
            edge3(1:end-yshift,:) = edge2(yshift+1:end,:);
        elseif (yshift < 0)
            img3 = zeros(size(img2,1),size(img2,2));
            img3(-yshift+1:end,:) = img2(1:end+yshift,:);
            img_2 = zeros(size(img2,1),size(img2,2));
            img_2(-yshift+1:end,:) = img_1(1:end+yshift,:);
            edge3 = zeros(size(img2,1),size(img2,2));
            edge3(-yshift+1:end,:) = edge2(1:end+yshift,:);
        else
            img3 = img2;
            img_2 = img_1;
            edge3 = edge2;
        end
        
        
      
        % Plot 
        set(h, 'AlphaData', transparency_green.*img3);
        % set(hhh, 'AlphaData', transparency.*edge2);
        set(hh, 'AlphaData', transparency_red.*I_REF); 
        set(min_h, 'AlphaData', overlay_bkg.*img_2);
        
        
        % plot_xy = plot(CELLS_X_ROT,CELLS_Y_ROT,'.g','Markersize',5);
        plot_xy = plot(CELLS_X_XSH,CELLS_Y_YSH,'.g','Markersize',5);
    end
    
    

    function slider_callback_x(hObj, eventdata)
        img = I;
        img2 = I;
        img3 = I;
        img_0 = MIN;
        img_1 = MIN;
        img_2 = MIN;
        edge1 = EDGES;
        edge2 = EDGES;
        edge2 = EDGES;
        reset(plot_xy);
        
        
%         rot = round(get(hObj,'Value'));   
%         set(hTxt_rot, 'String',strcat('Rotation: ',num2str(rot,'%.02f'))); 
        xshift = round(get(hObj,'Value'));
        set(hTxt_x, 'String', strcat('X shift: ',num2str(xshift,'%.02f')));
%         yshift = round(get(hObj,'Value'));
%         set(hTxt_y, 'String', strcat('Y shift: ',num2str(yshift,'%.02f')));
%         
%         transparency_green = get(hObj,'Value');
%         set(hTxt_transpg, 'String',strcat('Green: ',num2str(transparency_green,'%.02f')));
%         transparency_red = get(hObj,'Value');
%         set(hTxt_transpr, 'String',strcat('Red: ',num2str(transparency_red,'%.02f')));
%         overlay_bkg = get(hObj,'Value');
%         set(hTxt_transpb, 'String',strcat('BCK: ',num2str(overlay_bkg,'%.02f')));
%         overlay_ref = get(hObj,'Value');
%         set(hTxt_transpref, 'String',strcat('BCK: ',num2str(overlay_ref,'%.02f')));

        
        % Rotation
    
        img = imrotate(I,rot,'bilinear','crop'); 
        img_0 = imrotate(MIN,rot,'bilinear','crop');
        edge1 = imrotate(EDGES,rot,'bilinear','crop'); 
        edge2 = edge1;
        img2 = img;
        set(h, 'AlphaData', transparency_green.*img); 
        % set(hhh, 'AlphaData', transparency.*edge1); 
        set(hTxt_rot, 'String',strcat('Rotation: ',num2str(rot,'%.02f'))); 
        set(min_h, 'AlphaData', overlay_bkg.*img_0);    
        rot_matrix = [cos(-rot*pi/180) -sin(-rot*pi/180); sin(-rot*pi/180) cos(-rot*pi/180)];
        for k = 1:length(CELLS_X)
            [coord] = [width/2; height/2] + rot_matrix*([CELLS_X(k); CELLS_Y(k)]-[width/2; height/2]);
            CELLS_X_ROT(k) = coord(1);
            CELLS_Y_ROT(k) = coord(2);
        end
        
        
        
        % X shift
        
        CELLS_X_XSH = CELLS_X_ROT - xshift;
        if (xshift > 0)
            img2 = zeros(size(img,1),size(img,2));
            img2(:,1:end-xshift) = img(:,xshift+1:end);
            img_1 = zeros(size(img,1),size(img,2));
            img_1(:,1:end-xshift) = img_0(:,xshift+1:end);
            edge2 = zeros(size(img,1),size(img,2));
            edge2(:,1:end-xshift) = edge1(:,xshift+1:end);
        elseif (xshift < 0)
            img2 = zeros(size(img,1),size(img,2));
            img2(:,-xshift+1:end) = img(:,1:end+xshift);
            img_1 = zeros(size(img,1),size(img,2));
            img_1(:,-xshift+1:end) = img_0(:,1:end+xshift);
            edge2 = zeros(size(img,1),size(img,2));
            edge2(:,-xshift+1:end) = edge1(:,1:end+xshift);
        else
            img2 = img;
            edge2 = edge1;
            img_1 = img_0;
        end
    
        
        
        
        % Y shift
        
        CELLS_Y_YSH = CELLS_Y_ROT - yshift;
        
        if (yshift > 0)
            img3 = zeros(size(img2,1),size(img2,2));
            img3(1:end-yshift,:) = img2(yshift+1:end,:);
            img_2 = zeros(size(img2,1),size(img2,2));
            img_2(1:end-yshift,:) = img_1(yshift+1:end,:);
            edge3 = zeros(size(img2,1),size(img2,2));
            edge3(1:end-yshift,:) = edge2(yshift+1:end,:);
        elseif (yshift < 0)
            img3 = zeros(size(img2,1),size(img2,2));
            img3(-yshift+1:end,:) = img2(1:end+yshift,:);
            img_2 = zeros(size(img2,1),size(img2,2));
            img_2(-yshift+1:end,:) = img_1(1:end+yshift,:);
            edge3 = zeros(size(img2,1),size(img2,2));
            edge3(-yshift+1:end,:) = edge2(1:end+yshift,:);
        else
            img3 = img2;
            img_2 = img_1;
            edge3 = edge2;
        end

        
      
        % Plot 
        set(h, 'AlphaData', transparency_green.*img3);
        % set(hhh, 'AlphaData', transparency.*edge2);
        set(hh, 'AlphaData', transparency_red.*I_REF); 
        set(min_h, 'AlphaData', overlay_bkg.*img_2);
        
        
        % plot_xy = plot(CELLS_X_ROT,CELLS_Y_ROT,'.g','Markersize',5);
        plot_xy = plot(CELLS_X_XSH,CELLS_Y_YSH,'.g','Markersize',5);
    end


    

    function slider_callback_y(hObj, eventdata)
        img = I;
        img2 = I;
        img3 = I;
        img_0 = MIN;
        img_1 = MIN;
        img_2 = MIN;
        edge1 = EDGES;
        edge2 = EDGES;
        edge2 = EDGES;
        reset(plot_xy);
        
        
%         rot = round(get(hObj,'Value'));   
%         set(hTxt_rot, 'String',strcat('Rotation: ',num2str(rot,'%.02f'))); 
%         xshift = round(get(hObj,'Value'));
%         set(hTxt_x, 'String', strcat('X shift: ',num2str(xshift,'%.02f')));
        yshift = round(get(hObj,'Value'));
        set(hTxt_y, 'String', strcat('Y shift: ',num2str(yshift,'%.02f')));
%         
%         transparency_green = get(hObj,'Value');
%         set(hTxt_transpg, 'String',strcat('Green: ',num2str(transparency_green,'%.02f')));
%         transparency_red = get(hObj,'Value');
%         set(hTxt_transpr, 'String',strcat('Red: ',num2str(transparency_red,'%.02f')));
%         overlay_bkg = get(hObj,'Value');
%         set(hTxt_transpb, 'String',strcat('BCK: ',num2str(overlay_bkg,'%.02f')));
%         overlay_ref = get(hObj,'Value');
%         set(hTxt_transpref, 'String',strcat('BCK: ',num2str(overlay_ref,'%.02f')));

        
        % Rotation
    
        img = imrotate(I,rot,'bilinear','crop'); 
        img_0 = imrotate(MIN,rot,'bilinear','crop');
        edge1 = imrotate(EDGES,rot,'bilinear','crop'); 
        edge2 = edge1;
        img2 = img;
        set(h, 'AlphaData', transparency_green.*img); 
        % set(hhh, 'AlphaData', transparency.*edge1); 
        set(hTxt_rot, 'String',strcat('Rotation: ',num2str(rot,'%.02f'))); 
        set(min_h, 'AlphaData', overlay_bkg.*img_0);    
        rot_matrix = [cos(-rot*pi/180) -sin(-rot*pi/180); sin(-rot*pi/180) cos(-rot*pi/180)];
        for k = 1:length(CELLS_X)
            [coord] = [width/2; height/2] + rot_matrix*([CELLS_X(k); CELLS_Y(k)]-[width/2; height/2]);
            CELLS_X_ROT(k) = coord(1);
            CELLS_Y_ROT(k) = coord(2);
        end
        
        
        
        % X shift
        
        CELLS_X_XSH = CELLS_X_ROT - xshift;
        if (xshift > 0)
            img2 = zeros(size(img,1),size(img,2));
            img2(:,1:end-xshift) = img(:,xshift+1:end);
            img_1 = zeros(size(img,1),size(img,2));
            img_1(:,1:end-xshift) = img_0(:,xshift+1:end);
            edge2 = zeros(size(img,1),size(img,2));
            edge2(:,1:end-xshift) = edge1(:,xshift+1:end);
        elseif (xshift < 0)
            img2 = zeros(size(img,1),size(img,2));
            img2(:,-xshift+1:end) = img(:,1:end+xshift);
            img_1 = zeros(size(img,1),size(img,2));
            img_1(:,-xshift+1:end) = img_0(:,1:end+xshift);
            edge2 = zeros(size(img,1),size(img,2));
            edge2(:,-xshift+1:end) = edge1(:,1:end+xshift);
        else
            img2 = img;
            edge2 = edge1;
            img_1 = img_0;
        end
    
        
        
        
        % Y shift
        
        CELLS_Y_YSH = CELLS_Y_ROT - yshift;
        
        if (yshift > 0)
            img3 = zeros(size(img2,1),size(img2,2));
            img3(1:end-yshift,:) = img2(yshift+1:end,:);
            img_2 = zeros(size(img2,1),size(img2,2));
            img_2(1:end-yshift,:) = img_1(yshift+1:end,:);
            edge3 = zeros(size(img2,1),size(img2,2));
            edge3(1:end-yshift,:) = edge2(yshift+1:end,:);
        elseif (yshift < 0)
            img3 = zeros(size(img2,1),size(img2,2));
            img3(-yshift+1:end,:) = img2(1:end+yshift,:);
            img_2 = zeros(size(img2,1),size(img2,2));
            img_2(-yshift+1:end,:) = img_1(1:end+yshift,:);
            edge3 = zeros(size(img2,1),size(img2,2));
            edge3(-yshift+1:end,:) = edge2(1:end+yshift,:);
        else
            img3 = img2;
            img_2 = img_1;
            edge3 = edge2;
        end

        
      
        % Plot 
        set(h, 'AlphaData', transparency_green.*img3);
        % set(hhh, 'AlphaData', transparency.*edge2);
        set(hh, 'AlphaData', transparency_red.*I_REF); 
        set(min_h, 'AlphaData', overlay_bkg.*img_2);
        
        
        % plot_xy = plot(CELLS_X_ROT,CELLS_Y_ROT,'.g','Markersize',5);
        plot_xy = plot(CELLS_X_XSH,CELLS_Y_YSH,'.g','Markersize',5);
    end












    function slider_callback_transpg(hObj, eventdata)
        img = I;
        img2 = I;
        img3 = I;
        img_0 = MIN;
        img_1 = MIN;
        img_2 = MIN;
        edge1 = EDGES;
        edge2 = EDGES;
        edge2 = EDGES;
        reset(plot_xy);
        
        
%         rot = round(get(hObj,'Value'));   
%         set(hTxt_rot, 'String',strcat('Rotation: ',num2str(rot,'%.02f'))); 
%         xshift = round(get(hObj,'Value'));
%         set(hTxt_x, 'String', strcat('X shift: ',num2str(xshift,'%.02f')));
%         yshift = round(get(hObj,'Value'));
%         set(hTxt_y, 'String', strcat('Y shift: ',num2str(yshift,'%.02f')));
%         
        transparency_green = get(hObj,'Value');
        set(hTxt_transpg, 'String',strcat('Green: ',num2str(transparency_green,'%.02f')));
%         transparency_red = get(hObj,'Value');
%         set(hTxt_transpr, 'String',strcat('Red: ',num2str(transparency_red,'%.02f')));
%         overlay_bkg = get(hObj,'Value');
%         set(hTxt_transpb, 'String',strcat('BCK: ',num2str(overlay_bkg,'%.02f')));
%         overlay_ref = get(hObj,'Value');
%         set(hTxt_transpref, 'String',strcat('BCK: ',num2str(overlay_ref,'%.02f')));

        
        % Rotation
    
        img = imrotate(I,rot,'bilinear','crop'); 
        img_0 = imrotate(MIN,rot,'bilinear','crop');
        edge1 = imrotate(EDGES,rot,'bilinear','crop'); 
        edge2 = edge1;
        img2 = img;
        set(h, 'AlphaData', transparency_green.*img); 
        % set(hhh, 'AlphaData', transparency.*edge1); 
        set(hTxt_rot, 'String',strcat('Rotation: ',num2str(rot,'%.02f'))); 
        set(min_h, 'AlphaData', overlay_bkg.*img_0);    
        rot_matrix = [cos(-rot*pi/180) -sin(-rot*pi/180); sin(-rot*pi/180) cos(-rot*pi/180)];
        for k = 1:length(CELLS_X)
            [coord] = [width/2; height/2] + rot_matrix*([CELLS_X(k); CELLS_Y(k)]-[width/2; height/2]);
            CELLS_X_ROT(k) = coord(1);
            CELLS_Y_ROT(k) = coord(2);
        end
        
        
        
        % X shift
        
        CELLS_X_XSH = CELLS_X_ROT - xshift;
        if (xshift > 0)
            img2 = zeros(size(img,1),size(img,2));
            img2(:,1:end-xshift) = img(:,xshift+1:end);
            img_1 = zeros(size(img,1),size(img,2));
            img_1(:,1:end-xshift) = img_0(:,xshift+1:end);
            edge2 = zeros(size(img,1),size(img,2));
            edge2(:,1:end-xshift) = edge1(:,xshift+1:end);
        elseif (xshift < 0)
            img2 = zeros(size(img,1),size(img,2));
            img2(:,-xshift+1:end) = img(:,1:end+xshift);
            img_1 = zeros(size(img,1),size(img,2));
            img_1(:,-xshift+1:end) = img_0(:,1:end+xshift);
            edge2 = zeros(size(img,1),size(img,2));
            edge2(:,-xshift+1:end) = edge1(:,1:end+xshift);
        else
            img2 = img;
            edge2 = edge1;
            img_1 = img_0;
        end
  
        
        
        
        % Y shift
        
        CELLS_Y_YSH = CELLS_Y_ROT - yshift;
        
        if (yshift > 0)
            img3 = zeros(size(img2,1),size(img2,2));
            img3(1:end-yshift,:) = img2(yshift+1:end,:);
            img_2 = zeros(size(img2,1),size(img2,2));
            img_2(1:end-yshift,:) = img_1(yshift+1:end,:);
            edge3 = zeros(size(img2,1),size(img2,2));
            edge3(1:end-yshift,:) = edge2(yshift+1:end,:);
        elseif (yshift < 0)
            img3 = zeros(size(img2,1),size(img2,2));
            img3(-yshift+1:end,:) = img2(1:end+yshift,:);
            img_2 = zeros(size(img2,1),size(img2,2));
            img_2(-yshift+1:end,:) = img_1(1:end+yshift,:);
            edge3 = zeros(size(img2,1),size(img2,2));
            edge3(-yshift+1:end,:) = edge2(1:end+yshift,:);
        else
            img3 = img2;
            img_2 = img_1;
            edge3 = edge2;
        end

        
      
        % Plot 
        set(h, 'AlphaData', transparency_green.*img3);
        % set(hhh, 'AlphaData', transparency.*edge2);
        set(hh, 'AlphaData', transparency_red.*I_REF); 
        set(min_h, 'AlphaData', overlay_bkg.*img_2);
        
        
        % plot_xy = plot(CELLS_X_ROT,CELLS_Y_ROT,'.g','Markersize',5);
        plot_xy = plot(CELLS_X_XSH,CELLS_Y_YSH,'.g','Markersize',5);

    end

    function slider_callback_transpr(hObj, eventdata)
        img = I;
        img2 = I;
        img3 = I;
        img_0 = MIN;
        img_1 = MIN;
        img_2 = MIN;
        edge1 = EDGES;
        edge2 = EDGES;
        edge2 = EDGES;
        reset(plot_xy);
        
        
%         rot = round(get(hObj,'Value'));   
%         set(hTxt_rot, 'String',strcat('Rotation: ',num2str(rot,'%.02f'))); 
%         xshift = round(get(hObj,'Value'));
%         set(hTxt_x, 'String', strcat('X shift: ',num2str(xshift,'%.02f')));
%         yshift = round(get(hObj,'Value'));
%         set(hTxt_y, 'String', strcat('Y shift: ',num2str(yshift,'%.02f')));
%         
%         transparency_green = get(hObj,'Value');
%         set(hTxt_transpg, 'String',strcat('Green: ',num2str(transparency_green,'%.02f')));
        transparency_red = get(hObj,'Value');
        set(hTxt_transpr, 'String',strcat('Red: ',num2str(transparency_red,'%.02f')));
%         overlay_bkg = get(hObj,'Value');
%         set(hTxt_transpb, 'String',strcat('BCK: ',num2str(overlay_bkg,'%.02f')));
%         overlay_ref = get(hObj,'Value');
%         set(hTxt_transpref, 'String',strcat('BCK: ',num2str(overlay_ref,'%.02f')));

        
        % Rotation
    
        img = imrotate(I,rot,'bilinear','crop'); 
        img_0 = imrotate(MIN,rot,'bilinear','crop');
        edge1 = imrotate(EDGES,rot,'bilinear','crop'); 
        edge2 = edge1;
        img2 = img;
        set(h, 'AlphaData', transparency_green.*img); 
        % set(hhh, 'AlphaData', transparency.*edge1); 
        set(hTxt_rot, 'String',strcat('Rotation: ',num2str(rot,'%.02f'))); 
        set(min_h, 'AlphaData', overlay_bkg.*img_0);    
        rot_matrix = [cos(-rot*pi/180) -sin(-rot*pi/180); sin(-rot*pi/180) cos(-rot*pi/180)];
        for k = 1:length(CELLS_X)
            [coord] = [width/2; height/2] + rot_matrix*([CELLS_X(k); CELLS_Y(k)]-[width/2; height/2]);
            CELLS_X_ROT(k) = coord(1);
            CELLS_Y_ROT(k) = coord(2);
        end
        
        
        
        % X shift
        
        CELLS_X_XSH = CELLS_X_ROT - xshift;
        if (xshift > 0)
            img2 = zeros(size(img,1),size(img,2));
            img2(:,1:end-xshift) = img(:,xshift+1:end);
            img_1 = zeros(size(img,1),size(img,2));
            img_1(:,1:end-xshift) = img_0(:,xshift+1:end);
            edge2 = zeros(size(img,1),size(img,2));
            edge2(:,1:end-xshift) = edge1(:,xshift+1:end);
        elseif (xshift < 0)
            img2 = zeros(size(img,1),size(img,2));
            img2(:,-xshift+1:end) = img(:,1:end+xshift);
            img_1 = zeros(size(img,1),size(img,2));
            img_1(:,-xshift+1:end) = img_0(:,1:end+xshift);
            edge2 = zeros(size(img,1),size(img,2));
            edge2(:,-xshift+1:end) = edge1(:,1:end+xshift);
        else
            img2 = img;
            edge2 = edge1;
            img_1 = img_0;
        end
 
        
        
        
        % Y shift
        
        CELLS_Y_YSH = CELLS_Y_ROT - yshift;
        
        if (yshift > 0)
            img3 = zeros(size(img2,1),size(img2,2));
            img3(1:end-yshift,:) = img2(yshift+1:end,:);
            img_2 = zeros(size(img2,1),size(img2,2));
            img_2(1:end-yshift,:) = img_1(yshift+1:end,:);
            edge3 = zeros(size(img2,1),size(img2,2));
            edge3(1:end-yshift,:) = edge2(yshift+1:end,:);
        elseif (yshift < 0)
            img3 = zeros(size(img2,1),size(img2,2));
            img3(-yshift+1:end,:) = img2(1:end+yshift,:);
            img_2 = zeros(size(img2,1),size(img2,2));
            img_2(-yshift+1:end,:) = img_1(1:end+yshift,:);
            edge3 = zeros(size(img2,1),size(img2,2));
            edge3(-yshift+1:end,:) = edge2(1:end+yshift,:);
        else
            img3 = img2;
            img_2 = img_1;
            edge3 = edge2;
        end

        
      
        % Plot 
        set(h, 'AlphaData', transparency_green.*img3);
        % set(hhh, 'AlphaData', transparency.*edge2);
        set(hh, 'AlphaData', transparency_red.*I_REF); 
        set(min_h, 'AlphaData', overlay_bkg.*img_2);
        
        
        % plot_xy = plot(CELLS_X_ROT,CELLS_Y_ROT,'.g','Markersize',5);
        plot_xy = plot(CELLS_X_XSH,CELLS_Y_YSH,'.g','Markersize',5);

    end

    function slider_callback_transpb(hObj, eventdata)
        img = I;
        img2 = I;
        img3 = I;
        img_0 = MIN;
        img_1 = MIN;
        img_2 = MIN;
        edge1 = EDGES;
        edge2 = EDGES;
        edge2 = EDGES;
        reset(plot_xy);
        
        
%         rot = round(get(hObj,'Value'));   
%         set(hTxt_rot, 'String',strcat('Rotation: ',num2str(rot,'%.02f'))); 
%         xshift = round(get(hObj,'Value'));
%         set(hTxt_x, 'String', strcat('X shift: ',num2str(xshift,'%.02f')));
%         yshift = round(get(hObj,'Value'));
%         set(hTxt_y, 'String', strcat('Y shift: ',num2str(yshift,'%.02f')));
%         
%         transparency_green = get(hObj,'Value');
%         set(hTxt_transpg, 'String',strcat('Green: ',num2str(transparency_green,'%.02f')));
%         transparency_red = get(hObj,'Value');
%         set(hTxt_transpr, 'String',strcat('Red: ',num2str(transparency_red,'%.02f')));
        overlay_bkg = get(hObj,'Value');
        set(hTxt_transpb, 'String',strcat('BCK: ',num2str(overlay_bkg,'%.02f')));
%         overlay_ref = get(hObj,'Value');
%         set(hTxt_transpref, 'String',strcat('BCK: ',num2str(overlay_ref,'%.02f')));

        
        % Rotation
    
        img = imrotate(I,rot,'bilinear','crop'); 
        img_0 = imrotate(MIN,rot,'bilinear','crop');
        edge1 = imrotate(EDGES,rot,'bilinear','crop'); 
        edge2 = edge1;
        img2 = img;
        set(h, 'AlphaData', transparency_green.*img); 
        % set(hhh, 'AlphaData', transparency.*edge1); 
        set(hTxt_rot, 'String',strcat('Rotation: ',num2str(rot,'%.02f'))); 
        set(min_h, 'AlphaData', overlay_bkg.*img_0);    
        rot_matrix = [cos(-rot*pi/180) -sin(-rot*pi/180); sin(-rot*pi/180) cos(-rot*pi/180)];
        for k = 1:length(CELLS_X)
            [coord] = [width/2; height/2] + rot_matrix*([CELLS_X(k); CELLS_Y(k)]-[width/2; height/2]);
            CELLS_X_ROT(k) = coord(1);
            CELLS_Y_ROT(k) = coord(2);
        end
        
        
        
        % X shift
        
        CELLS_X_XSH = CELLS_X_ROT - xshift;
        if (xshift > 0)
            img2 = zeros(size(img,1),size(img,2));
            img2(:,1:end-xshift) = img(:,xshift+1:end);
            img_1 = zeros(size(img,1),size(img,2));
            img_1(:,1:end-xshift) = img_0(:,xshift+1:end);
            edge2 = zeros(size(img,1),size(img,2));
            edge2(:,1:end-xshift) = edge1(:,xshift+1:end);
        elseif (xshift < 0)
            img2 = zeros(size(img,1),size(img,2));
            img2(:,-xshift+1:end) = img(:,1:end+xshift);
            img_1 = zeros(size(img,1),size(img,2));
            img_1(:,-xshift+1:end) = img_0(:,1:end+xshift);
            edge2 = zeros(size(img,1),size(img,2));
            edge2(:,-xshift+1:end) = edge1(:,1:end+xshift);
        else
            img2 = img;
            edge2 = edge1;
            img_1 = img_0;
        end
   
        
        
        
        % Y shift
        
        CELLS_Y_YSH = CELLS_Y_ROT - yshift;
        
        if (yshift > 0)
            img3 = zeros(size(img2,1),size(img2,2));
            img3(1:end-yshift,:) = img2(yshift+1:end,:);
            img_2 = zeros(size(img2,1),size(img2,2));
            img_2(1:end-yshift,:) = img_1(yshift+1:end,:);
            edge3 = zeros(size(img2,1),size(img2,2));
            edge3(1:end-yshift,:) = edge2(yshift+1:end,:);
        elseif (yshift < 0)
            img3 = zeros(size(img2,1),size(img2,2));
            img3(-yshift+1:end,:) = img2(1:end+yshift,:);
            img_2 = zeros(size(img2,1),size(img2,2));
            img_2(-yshift+1:end,:) = img_1(1:end+yshift,:);
            edge3 = zeros(size(img2,1),size(img2,2));
            edge3(-yshift+1:end,:) = edge2(1:end+yshift,:);
        else
            img3 = img2;
            img_2 = img_1;
            edge3 = edge2;
        end

        
      
        % Plot 
        set(h, 'AlphaData', transparency_green.*img3);
        % set(hhh, 'AlphaData', transparency.*edge2);
        set(hh, 'AlphaData', transparency_red.*I_REF); 
        set(min_h, 'AlphaData', overlay_bkg.*img_2);
        
        
        % plot_xy = plot(CELLS_X_ROT,CELLS_Y_ROT,'.g','Markersize',5);
        plot_xy = plot(CELLS_X_XSH,CELLS_Y_YSH,'.g','Markersize',5);
    end



end