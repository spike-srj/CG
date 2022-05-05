void rst::rasterizer::rasterize_triangle(const Triangle& t, const std::array<Eigen::Vector3f, 3>& view_pos) 
{
    auto v = t.toVector4();

    std::vector<float> x_arry{ v[0].x(), v[1].x(), v[2].x() };
    std::vector<float> y_arry{ v[0].y(), v[1].y(), v[2].y() };
    std::sort(x_arry.begin(), x_arry.end());
    std::sort(y_arry.begin(), y_arry.end());
    int x_min = floor(x_arry[0]), x_max = ceil(x_arry[2]),
        y_min = floor(y_arry[0]), y_max = ceil(y_arry[2]);

    for (int x = x_min; x < x_max; x++)
    {
        for (int y = y_min; y < y_max; y++) {

            Eigen::Vector2i point(x, y); //这个注意是2i，之前是3f，不改过来会有奇怪的错误、

                if (insideTriangle(x + 0.5f, y + 0.5f, t.v)) {
                    auto [alpha, beta, gamma] = computeBarycentric2D(x + 0.5, y + 0.5, t.v);
                    float w_reciprocal = 1.0 / (alpha / v[0].w() + beta / v[1].w() + gamma / v[2].w());
                    float z_interpolated = alpha * v[0].z() / v[0].w() + beta * v[1].z() / v[1].w() + gamma * v[2].z() / v[2].w();
                    z_interpolated *= w_reciprocal;

                    auto normal_interpolated = alpha * t.normal[0] + beta * t.normal[1] + gamma * t.normal[2];
                    auto color_interpolated = alpha * t.color[0] + beta * t.color[1] + gamma * t.color[2];
                    auto textureCoord_interpolated = alpha * t.tex_coords[0] + beta * t.tex_coords[1] + gamma * t.tex_coords[2];
                    auto shadingcoords_interpolated = alpha * view_pos[0] + beta * view_pos[1] + gamma * view_pos[2];

                    if (z_interpolated < depth_buf[get_index(x, y)]) {
                        fragment_shader_payload payload(color_interpolated,normal_interpolated.normalized(), textureCoord_interpolated, texture ? &*texture : nullptr);
                        payload.view_pos = shadingcoords_interpolated;
                        auto pixel_color = fragment_shader(payload);
                        set_pixel(point, pixel_color);
                        depth_buf[get_index(x, y)] = z_interpolated;
                   }
                }
        }
    }
}


















void rst::rasterizer::rasterize_triangle(const Triangle& t, const std::array<Eigen::Vector3f, 3>& view_pos) 
{
    // TODO: From your HW3, get the triangle rasterization code.
    auto v = t.toVector4();
    float x_min = std::min(v[0].x(),std::min(v[1].x(),v[2].x()));
    float y_min = std::min(v[0].y(),std::min(v[1].y(),v[2].y()));
    float x_max = std::max(v[0].x(),std::max(v[1].x(),v[2].x()));
    float y_max = std::max(v[0].y(),std::max(v[1].y(),v[2].y()));
    // TODO : Find out the bounding box of current triangle.
    // iterate through the pixel and find if the current pixel is inside the triangle
    for (int x=x_min; x<x_max; x++)
    {
        for (int y = y_min; y < y_max; y++)
        {
            if (insideTriangle(x+0.5f, y+0.5f, t.v))
            {
                float alpha, beta, gamma;
                auto tup = computeBarycentric2D(x, y, t.v);
                std::tie(alpha,beta,gamma) = tup;
                float w_reciprocal = 1.0/(alpha / v[0].w() + beta / v[1].w() + gamma / v[2].w());
                float z_interpolated = alpha * v[0].z() / v[0].w() + beta * v[1].z() / v[1].w() + gamma * v[2].z() / v[2].w();
                z_interpolated *= w_reciprocal;
                // TODO : set the current pixel (use the set_pixel function) to the color of the triangle (use getColor function) if it should be painted.
                int index = get_index(x,y);
                auto interpolated_color = alpha*t.color[0]+beta*t.color[1]+gamma*t.color[2];
                auto interpolated_normal = alpha*t.normal[0]+beta*t.normal[1]+gamma*t.normal[2];
                auto interpolated_texcoords = alpha*t.tex_coords[0]+beta*t.tex_coords[1]+gamma*t.tex_coords[2];
                auto interpolated_shadingcoords = alpha*view_pos[0]+beta*view_pos[1]+gamma*view_pos[2];
                if (z_interpolated < depth_buf[index])//
                {
                    
                    fragment_shader_payload payload(interpolated_color,interpolated_normal.normalized(),interpolated_texcoords,texture ?&*texture : nullptr);
                    payload.view_pos = interpolated_shadingcoords;
                    auto pixel_color = fragment_shader(payload);
                    Eigen::Vector2i p;
                    p << x,y;
                    set_pixel(p,pixel_color);
                    depth_buf[index] = z_interpolated;
                }
            }
        }
        
    }
    // TODO: Inside your rasterization loop:
    //    * v[i].w() is the vertex view space depth value z.
    //    * Z is interpolated view space depth for the current pixel
    //    * zp is depth between zNear and zFar, used for z-buffer

    // float Z = 1.0 / (alpha / v[0].w() + beta / v[1].w() + gamma / v[2].w());
    // float zp = alpha * v[0].z() / v[0].w() + beta * v[1].z() / v[1].w() + gamma * v[2].z() / v[2].w();
    // zp *= Z;

    // TODO: Interpolate the attributes:
    // auto interpolated_color
    // auto interpolated_normal
    // auto interpolated_texcoords
    // auto interpolated_shadingcoords

    // Use: fragment_shader_payload payload( interpolated_color, interpolated_normal.normalized(), interpolated_texcoords, texture ? &*texture : nullptr);
    // Use: payload.view_pos = interpolated_shadingcoords;
    // Use: Instead of passing the triangle's color directly to the frame buffer, pass the color to the shaders first to get the final color;
    // Use: auto pixel_color = fragment_shader(payload);

 
}