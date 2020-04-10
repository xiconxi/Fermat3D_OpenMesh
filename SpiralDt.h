//
// Created by hotpot on 2020/4/5.
//

#ifndef VISMAGNUM_SPIRALDT_H
#define VISMAGNUM_SPIRALDT_H

#include <Eigen/Core>
#include <vector>

#include <ostream>


struct Handle{
    int i;
    bool is_valid(){return i != -1;}
};

struct LinkNode{
    int prev, next;
    LinkNode(const LinkNode&)= delete;
    LinkNode(LinkNode&&) = default;
};


struct VLinkNode: LinkNode{
    size_t iso_i;
};

struct ILinkNode: LinkNode{
    int parent, child;
    int iso_level, v_i;
};

struct IterLine{
    struct VHandle: Handle{};
    struct IHandle: Handle{};
    struct RegionNode {
        IHandle begin, end;
        int parent;
    };

    VHandle next(VHandle vh) const{ return VHandle{Vs[vh.i].next}; }
    VHandle prev(VHandle vh) const{ return VHandle{Vs[vh.i].prev}; }
    IHandle next(IHandle ih) const{ return IHandle{Is[ih.i].next}; }
    IHandle prev(IHandle ih) const{ return IHandle{Is[ih.i].prev}; }
    IHandle parent(IHandle ih) const{ return IHandle{Is[ih.i].parent}; }
    IHandle child(IHandle ih) const{ return IHandle{Is[ih.i].child}; }
    VHandle v_handle(IHandle ih) const{ return VHandle{int(Is[ih.i].v_i)}; }
    IHandle i_handle(VHandle vh) const{ return IHandle{int(Vs[vh.i].iso_i)}; }
    VLinkNode& link(VHandle vh) { return Vs[vh.i]; }
    ILinkNode& iso(IHandle ih) { return Is[ih.i]; }
    bool is_isolated (VHandle vh) const{ return Vs[vh.i].next == -1 and Vs[vh.i].prev == -1; }
    double length(IHandle ih);

    VHandle add_isoline(const std::vector<Eigen::RowVector3d>& iso_v);
    
    IHandle add_hierarchy(VHandle vh, size_t iso_level);

    Eigen::RowVector3d& p(VHandle vh)  { return  Ps[vh.i]; }

    double hausdorff(IHandle lhs, IHandle rhs);
    VHandle hausdorff(VHandle vh, IHandle ih);
    VHandle hausdorff(VHandle vh, VHandle spiral_begin);

    bool after_test(VHandle vh, VHandle next_vh);

    void after_delete(VHandle s_vh, VHandle e_vh);
    VHandle before(VHandle vh, double length);
    VHandle after(VHandle vh, double length);
    void after_reverse(VHandle s_vh, VHandle e_vh);

    std::vector<RegionNode> region_extract();
    std::tuple<VHandle,VHandle>  make_fermat_spiral(std::vector<RegionNode>& regions);

    void print_fermat_spiral(std::ostream& stream, VHandle vh);

    std::tuple<VHandle,VHandle> spiral_from_maximum(IHandle s_ih, IHandle e_ih);
 
private:

    std::vector<Eigen::RowVector3d> Ps;
    std::vector<VLinkNode> Vs;
    std::vector<ILinkNode> Is;

    friend std::ostream& operator << (std::ostream& stream, const IterLine&);
};

std::ostream& operator << (std::ostream& stream, const IterLine&);


#endif //VISMAGNUM_SPIRALDT_H
