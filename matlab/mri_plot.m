function plot_tracts_with_mri()

    %  Click on any the three MRI images and press "ENTER" to change slices at new [x,y,z] positions.
    %  Press "ESC" and press "ENTER" to exit.
    %  by Binlin Wu -- CCNY
    %  bwu@sci.ccny.cuny.edu
    %  09/14/2010 

    load mri.mat;
    D1=double(squeeze(D));
    DIM = size(D1);
    [X,Y,Z]=meshgrid(1:DIM(2),1:DIM(1),1:DIM(3));
    h1=subplot(2,2,1);imagesc(D1(:,:,round(DIM(3)/2)),[min(D1(:)) max(D1(:))]);colormap(gray);title('axial');colorbar;
    xlabel('x');ylabel('y')
    h2=subplot(2,2,2);imagesc(squeeze(D1(:,round(DIM(2)/2),:)),[min(D(:)) max(D(:))]);colormap(gray);title('sagittal');colorbar;
    xlabel('z');ylabel('y')
    h3=subplot(2,2,3);imagesc(squeeze(D1(round(DIM(1)/2),:,:)),[min(D(:)) max(D(:))]);colormap(gray);title('coronal');colorbar;
    xlabel('z');ylabel('x')
    subplot(2,2,4);slice(X,Y,Z,D1,64,64,14);colormap(gray);shading flat;title('3D Slices')
    xlabel('x');ylabel('y');zlabel('z');


    x=round(DIM(2)/2);y=round(DIM(1)/2);z=round(DIM(3)/2);
    button = 0;

    while(1)

        try
            [A,B,button]=ginput
        catch
            return
        end
        if length(A)==0
            A=14;B=64;button=0;
        end

        A=A(end);
        B=B(end);
        button=button(end);

        A=ceil(A-0.5);
        B=ceil(B-0.5);

        if button==27
            break;
        end
        if gca==h1
            x=A;
            y=B;
            if x<=0 || x>DIM(2) || y<=0 || y>DIM(1)
                continue
            end
            axes(h2);imagesc(squeeze(D1(:,x,:)),[min(D(:)) max(D1(:))]);colormap(gray);title('sagittal');colorbar;
            xlabel('z');ylabel('y')
            axes(h3);imagesc(squeeze(D1(y,:,:)),[min(D(:)) max(D1(:))]);colormap(gray);title('coronal');colorbar;
            xlabel('z');ylabel('x')
            subplot(2,2,4);slice(X,Y,Z,D1,x,y,z);colormap(gray);shading flat;title('3D Slices')
            xlabel('x');ylabel('y');zlabel('z');
        elseif gca==h2
            z=A;
            y=B;
            if z<=0 || z>DIM(3) || y<=0 || y>DIM(1)
                continue
            end        
            axes(h1);imagesc(D1(:,:,z),[min(D1(:)) max(D1(:))]);colormap(gray);title('axial');colorbar;
            xlabel('x');ylabel('y')
            axes(h3);imagesc(squeeze(D1(y,:,:)),[min(D1(:)) max(D1(:))]);colormap(gray);title('coronal');colorbar;
            xlabel('z');ylabel('x')
            subplot(2,2,4);slice(X,Y,Z,D1,x,y,z);colormap(gray);shading flat;title('3D Slices')
            xlabel('x');ylabel('y');zlabel('z');        
        elseif gca==h3
            z=A;
            x=B;
            if x<=0 || x>DIM(2) || z<=0 || z>DIM(3)
                continue
            end        
            axes(h1);imagesc(D1(:,:,z),[min(D1(:)) max(D1(:))]);colormap(gray);title('axial');colorbar;
            xlabel('x');ylabel('y')
            axes(h2);imagesc(squeeze(D1(:,x,:)),[min(D1(:)) max(D1(:))]);colormap(gray);title('sagittal');colorbar;
            xlabel('z');ylabel('y')
            subplot(2,2,4);slice(X,Y,Z,D1,x,y,z);colormap(gray);shading flat;title('3D Slices')
            xlabel('x');ylabel('y');zlabel('z');        
        end
    end