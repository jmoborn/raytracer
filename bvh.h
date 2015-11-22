#ifndef __BVH_H__
#define __BVH_H__

#include <cstddef>
#include <iostream>
#include <vector>
#include <algorithm>
#include "vec4.h"

// template class must have bounding box variables (cuz i'm lazy)
template <class T>
class bvh_node
{
public:
    bvh_node()
    {
        this->left = NULL;
        this->right = NULL;
    };
    ~bvh_node(){};

    bvh_node<T> *left;
    bvh_node<T> *right;
    std::vector<T*> prims;

    double min[3];
    double max[3];
    // double xmin;
    // double xmax;
    // double ymin;
    // double ymax;
    // double zmin;
    // double zmax;
};

template <class T>
class bvh
{
public:
    bvh()
    {
        root = NULL;
    }
    ~bvh(){};
    void build(std::vector<T*> scene, double xmin, double xmax, double ymin, double ymax, double zmin, double zmax)
    {
        root = new bvh_node<T>;
        root->min[0] = xmin;
        root->min[1] = ymin;
        root->min[2] = zmin;
        root->max[0] = xmax;
        root->max[1] = ymax;
        root->max[2] = zmax;
        root->prims = scene;
        build_bvh(root, 1);
    };

    bvh_node<T> *root;
private:
    void build_bvh(bvh_node<T> *node, int level)
    {
        if(node->prims.size() <=10 || level >= 20)
        {
            // std::cout << "bvh node: " << node->prims.size() << " prims" << std::endl;
            // std::cout << "    level: " << level << std::endl;
            return;
        }
        double size[3];
        size[0] = node->max[0] - node->min[0];
        size[1] = node->max[1] - node->min[1];
        size[2] = node->max[2] - node->min[2];
        int largest_dim = -1;
        if (size[0] >= size[1] && size[0] >= size[2])
        {
            largest_dim = 0;
        }
        else if(size[1] >= size[0] && size[1] >= size[2])
        {
            largest_dim = 1;
        }
        else
        {
            largest_dim = 2;
        }
        // if(level==2||level==1) std::cout << "largest_dim: " << largest_dim << std::endl;

        bvh_node<T> *left_node = new bvh_node<T>;
        bvh_node<T> *right_node = new bvh_node<T>;
        for(int i=0; i<3; i++)
        {
            left_node->min[i] = node->min[i];
            right_node->max[i] = node->max[i];
            if(i==largest_dim)
            {
                left_node->max[i] = right_node->min[i] = node->max[i] - size[i]/2.0;
                // right_node->min[i] = node->min[i] + size[i]/2.0;
            }
            else
            {
                left_node->max[i] = node->max[i];
                right_node->min[i] = node->min[i];
            }
        }
//         std::cout << "parent bbox min: " << node->min[0] << " " << node->min[1] << " " << node->min[2] << std::endl;
//         std::cout << "parent bbox max: " << node->max[0] << " " << node->max[1] << " " << node->max[2] << std::endl;

        // std::cout << "left bbox min: " << left_node->min[0] << " " << left_node->min[1] << " " << left_node->min[2] << std::endl;
        // std::cout << "left bbox max: " << left_node->max[0] << " " << left_node->max[1] << " " << left_node->max[2] << std::endl;

//         std::cout << "right bbox min: " << right_node->min[0] << " " << right_node->min[1] << " " << right_node->min[2] << std::endl;
//         std::cout << "right bbox max: " << right_node->max[0] << " " << right_node->max[1] << " " << right_node->max[2] << std::endl;

// std::cout << "prims: " << node->prims.size() << std::endl;
        for(int i=0; i<node->prims.size(); i++)
        {
// std::cout << "prim min: " << node->prims[i]->xmin << " " << node->prims[i]->ymin << " " << node->prims[i]->zmin << std::endl;
// std::cout << "prim max: " << node->prims[i]->xmax << " " << node->prims[i]->ymax << " " << node->prims[i]->zmax << std::endl;
            if((node->prims[i]->xmin >= left_node->min[0] && node->prims[i]->xmin <= left_node->max[0]) ||
               (node->prims[i]->xmax >= left_node->min[0] && node->prims[i]->xmax <= left_node->max[0]) ||
               (node->prims[i]->xmin <= left_node->min[0] && node->prims[i]->xmax >= left_node->max[0]))
            {
                if((node->prims[i]->ymin >= left_node->min[1] && node->prims[i]->ymin <= left_node->max[1]) ||
                   (node->prims[i]->ymax >= left_node->min[1] && node->prims[i]->ymax <= left_node->max[1]) ||
                   (node->prims[i]->ymin <= left_node->min[1] && node->prims[i]->ymax >= left_node->max[1]))
                {
                    if((node->prims[i]->zmin >= left_node->min[2] && node->prims[i]->zmin <= left_node->max[2]) ||
                       (node->prims[i]->zmax >= left_node->min[2] && node->prims[i]->zmax <= left_node->max[2]) ||
                       (node->prims[i]->zmin <= left_node->min[2] && node->prims[i]->zmax >= left_node->max[2]))
                    {
                        left_node->prims.push_back(node->prims[i]);
                    }
                }
            }
            if((node->prims[i]->xmin >= right_node->min[0] && node->prims[i]->xmin <= right_node->max[0]) ||
               (node->prims[i]->xmax >= right_node->min[0] && node->prims[i]->xmax <= right_node->max[0]) ||
               (node->prims[i]->xmin <= right_node->min[0] && node->prims[i]->xmax >= right_node->max[0]))
            {
                if((node->prims[i]->ymin >= right_node->min[1] && node->prims[i]->ymin <= right_node->max[1]) ||
                   (node->prims[i]->ymax >= right_node->min[1] && node->prims[i]->ymax <= right_node->max[1]) ||
                   (node->prims[i]->ymin <= right_node->min[1] && node->prims[i]->ymax >= right_node->max[1]))
                {
                    if((node->prims[i]->zmin >= right_node->min[2] && node->prims[i]->zmin <= right_node->max[2]) ||
                       (node->prims[i]->zmax >= right_node->min[2] && node->prims[i]->zmax <= right_node->max[2]) ||
                       (node->prims[i]->zmin <= right_node->min[2] && node->prims[i]->zmax >= right_node->max[2]))
                    {
                        right_node->prims.push_back(node->prims[i]);
                    }
                }
            }
        }
        node->left = left_node;
        node->right = right_node;
        if((left_node->prims.size() + right_node->prims.size())<node->prims.size())
        {
            std::cout << "error: missing polygons!!!" << std::endl;
        }
        build_bvh(left_node, level+1);
        build_bvh(right_node, level+1);       
    }
};

#endif