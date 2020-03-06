% plot_z_surface. 
% By David Dorran (david.dorran@dit.ie)
%
%
 
function plot_z_surface(ax, pole_positions_orig, zero_positions_orig, varargin)
surface_limit = 1.5;
min_val = -surface_limit;
max_val = surface_limit;

CameraPos=[0 0 519.6152];
CameraUpVec=[0 0 1];
if(nargin >3 )
   CameraPos = varargin{1};
   CameraUpVec = varargin{2};
end
surface_display_opts = 0;
if(nargin == 5)
    surface_display_opts = varargin{3};
end
 
%check if any of the poles or zeros fall outside the maximum values
%specified
if length(find(abs(real([pole_positions_orig zero_positions_orig])) > max_val)) || length(find(abs(imag([pole_positions_orig zero_positions_orig])) > max_val))
    error(['the real part of the complex numbers indicating the the positions of the poles and zeros must be less than ' num2str(max_val) ' and greater than -' num2str(max_val) '. Also, the imaginary part of the complex numbers indicating the the positions of the poles and zeros must be less than ' num2str(max_val) 'j and greater than -' num2str(max_val) 'j.']);
 
end
 
%set up the s_plane grid from -5j to 5j along the imaginary axis and -5 to
%5 along the real axis. A resolution of 0.02 provides a visually pleasing
%plot
surface_resolution = 0.02;
 
a_vals = min_val:surface_resolution:max_val;
b_vals = a_vals;
grid_len = length(a_vals);
%create a grid of values along the z surface indicating position on
%the z-domain
ones_vector = ones(1, grid_len);
z_grid = ones_vector'*a_vals + (fliplr(b_vals*j))'*ones_vector;
 
%round all the poles and zeros to 1 decimal place because of resolution limitation
%Also, add .01+0.01j to each pole and zero position so that they are positioned in the
%center of a grid square of the z-surface created above. This ensures that
%the all the poles and zeros appear the same on the z-surface given the
%resolution constraints
pole_positions = round(pole_positions_orig*10)/10 + surface_resolution/2+(surface_resolution/2)*j;
zero_positions = round(zero_positions_orig*10)/10 +surface_resolution/2 +(surface_resolution/2)*j;
 
z_surface = zeros(length(a_vals)); % initialise the s_surface
% Since zeros are numerator terms obtain an log s_plane surface for each zero and add each individual plane
%  since a logarithmic addition is equivalent to a multiply i.e.
% log(a.b) = log(a)+log(b)
 
for k = 1 : length(zero_positions)
    z_surface = z_surface + 20*log10(abs(z_grid - zero_positions(k)));
%    transfer_function_numerator = [transfer_function_numerator '(s-' num2str(zero_positions_orig(k)),')'];
end
 
% Since poles are denominator terms obtain an log s_plane surface for each pole and subtract each individual plane
%  since a logarithmic subtraction is equivalent to a divide i.e.
% log(a/b) = log(a)-log(b)
transfer_function_denominator = [];
for k = 1 :length(pole_positions)
    z_surface = z_surface - 20*log10(abs(z_grid - pole_positions(k) ));
 
    transfer_function_denominator = [transfer_function_denominator '(z-' num2str(pole_positions_orig(k)),')'];
end
% 
% for k = 1 :length(pole_positions)
%     round(imag(pole_positions(k))/surface_resolution)
%     z_surface(round(-imag(pole_positions(k))/surface_resolution)+surface_limit/surface_resolution, round(real(pole_positions(k))/surface_resolution)+surface_limit/surface_resolution) = 40;
% end
% for k = 1 :length(zero_positions)
%     z_surface(round(imag(zero_positions(k))/surface_resolution)+surface_limit/surface_resolution, round(real(zero_positions(k))/surface_resolution)+surface_limit/surface_resolution) = -40;
% end
 
%half_z_surface = s_surface;
%half_s_surface(:,1:round(length(a_vals)/2)) = 20*log10(0.001);
 
 
 
% transfer_function = num2str(gain);
% if(length(transfer_function_numerator))
%     if(gain ~= 1)
%         transfer_function = [num2str(gain) '(' transfer_function_numerator ')'];
%     else
%         transfer_function = transfer_function_numerator;
%     end
% end
 
% if(length(transfer_function_denominator))
%     transfer_function = [ transfer_function '/(' transfer_function_denominator ')'];
% end
 
 
total_indices = 1:numel(z_grid);
temp_zero = 20*log10(abs(z_grid+.000001));
indices_outside_unit_circle = find(temp_zero > -0.2);
indices_inside_unit_circle = find(temp_zero < 0.2);
indices_on_unit_circle = intersect(indices_outside_unit_circle, indices_inside_unit_circle);
indices_not_on_unit_circle = setdiff( total_indices, indices_on_unit_circle );
positive_imaginary_indices = find(imag(z_grid) > 0 );
negative_inaginary_indices = find(imag(z_grid) < 0 );
positive_indices_outside_unit_circle = intersect(positive_imaginary_indices, indices_outside_unit_circle);
 
 
mask = ones(size(z_surface));
if(surface_display_opts == 1)
    mask_indices = indices_not_on_unit_circle;
     
else
    mask_indices = [];
end
mask(mask_indices) = NaN;
colors = z_surface;
colors(indices_on_unit_circle) = 40;
surf(a_vals,b_vals,z_surface.*mask, colors, 'EdgeColor','none')
set(gca,'YTick',[min_val:max_val])
set(gca,'YTickLabel',[min_val:max_val]')
set(gca,'XTick',[min_val:max_val])
set(gca,'XTickLabel',[min_val:max_val]')
campos(ax,CameraPos);
camup(ax,CameraUpVec);
%set(gca,'CameraPosition',CameraPos);
%set(gca,'CameraUpVector',CameraUpVec);
xlh = ylabel('Im');
ylh = xlabel('Re');
zlh = zlabel('Mag(db)');
axis([-1.1 1.1 -1.1 1.1 -30 30]) 
axis square
box on