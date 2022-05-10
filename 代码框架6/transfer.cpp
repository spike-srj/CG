    Vector3f invdir(1./ray.direction.x,1./ray.direction.y,1./ray.direction.z);
    std::array<int, 3> dirIsNeg;
    dirIsNeg[0] = ray.direction.x>0;
    dirIsNeg[1] = ray.direction.y>0;
    dirIsNeg[2] = ray.direction.z>0;
     if(!node->bounds.IntersectP(ray,invdir,dirIsNeg))
     {
         return intersect;
     }
    if(node->left == nullptr && node->right==nullptr)
    {
       return node->object->getIntersection(ray);
    }
    Intersection h1 = getIntersection(node->left,ray);
    Intersection h2 = getIntersection(node->right,ray);
    return h1.distance<h2.distance?h1:h2;
    return intersect;