function display_reference_sphere(sphere_radius, face_alpha)

  centre = [0; 0; 0];

  num_faces = round(sphere_radius * 100);

  [X, Y, Z] = sphere(num_faces);

  X = X * sphere_radius + centre(1);
  Y = Y * sphere_radius + centre(2);
  Z = Z * sphere_radius + centre(3);

  hold on;    

  h = surf(X, Y, Z);

  hold off;

  set(h,'facecolor', [0.5 0.5 0.5]);
  set(h, 'FaceAlpha', face_alpha);
  set(h,'edgecolor', 'none')


end