#pragma once

#include "typedef.h"
#include <list>
// Container class for "ff", the big field which contains the objects manipulated in the solvers.


class MeshContainer {
public:
        std::list<Mesh> F_;

        // Default constructor.
        MeshContainer(std::list<Mesh> InputMeshes)
        {
                F_ = InputMeshes;
        }

}; // End class MeshContainer.
