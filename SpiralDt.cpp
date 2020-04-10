//
// Created by hotpot on 2020/4/5.
//

#include "SpiralDt.h"

#include <iostream>

using VHandle = IterLine::VHandle;
using IHandle = IterLine::IHandle;
using RegionNode = IterLine::RegionNode;

std::ostream& operator << (std::ostream& stream, const IterLine& L) {
    for(int i = 0; i < L.Vs.size(); i++)
        stream << "v " << L.Ps[i][0] << ' ' << L.Ps[i][1]  <<' ' << L.Ps[i][2]  << std::endl;
    for(int i = 0; i < L.Vs.size(); i++) {
        if(L.is_isolated(IterLine::VHandle{i})) continue;
        if(L.Vs[i].next  == -1) continue;
        stream << "l " << i + 1 << ' ' << L.Vs[i].next + 1 << std::endl;
    }
    return  stream;
}

VHandle IterLine::add_isoline(const std::vector<Eigen::RowVector3d>& iso_v){
    int old_pts = Vs.size();
    Ps.insert(Ps.end(), iso_v.begin(), iso_v.end());

    for(int i = 0; i < iso_v.size() ; i++)
        Vs.push_back(VLinkNode{old_pts + i - 1, old_pts + i + 1, Is.size()});

    Vs[old_pts].prev = Vs.size()-1;
    Vs[Vs.size()-1].next = old_pts;
    return VHandle{int(Vs.size()-1)};
}

IHandle IterLine::add_hierarchy(VHandle vh, size_t iso_level) {
    Is.push_back(ILinkNode{-1, -1, -1, -1, int(iso_level), vh.i});
    double min_hausdorff = 1e10;
    int parent_i = -1;
    for(int i = 0; i < Is.size(); i++) {
        if(Is[i].iso_level + 1 != iso_level) continue;
        double _hausdorff = hausdorff(IHandle{i}, IHandle{int(Is.size()-1)});
        if( _hausdorff < min_hausdorff) {
            min_hausdorff = _hausdorff;
            parent_i = i;
        }
    }

    if (parent_i == -1)
        return IHandle{int(Is.size()-1)};

    if(Is[parent_i].child != -1) {
        Is.rbegin()->next = Is[parent_i].child;
        Is[Is[parent_i].child].prev = Is.size() - 1;
    }

    Is.rbegin()->parent = parent_i;
    Is[parent_i].child = Is.size() - 1;

    return IHandle{int(Is.size()-1)};
}

double IterLine::hausdorff(IHandle lh, IHandle rh) {
    assert(labs( Is[lh.i].iso_level-Is[rh.i].iso_level) == 1);
    VHandle lhs_s = v_handle(lh), rhs_s = v_handle(rh);
    double hausdorff = (p(lhs_s)-p(rhs_s)).norm();
    // trick implement, both lhs_s and rhs_s  are not included.
    for(VHandle i = next(lhs_s); i.i != lhs_s.i; i = next(i))
        for(VHandle j = next(rhs_s); j.i != rhs_s.i; j = next(j))
            hausdorff = std::min(hausdorff, (p(i)-p(j)).norm());

    for(VHandle j = next(rhs_s); j.i != rhs_s.i; j = next(j))
        hausdorff = std::min(hausdorff, (p(lhs_s)-p(j)).norm());
    for(VHandle i = next(lhs_s); i.i != lhs_s.i; i = next(i))
        hausdorff = std::min(hausdorff, (p(i)-p(rhs_s)).norm());

    return hausdorff;
}

VHandle IterLine::hausdorff(VHandle vh, IHandle ih) {
    VHandle t_s = v_handle(ih);
    double hausdorff = (p(vh)-p(t_s)).norm();
    VHandle result = t_s;

    for(VHandle j = next(t_s); j.i != t_s.i; j = next(j)) {
        double d = (p(vh) - p(j)).norm();
        if (d  < hausdorff) {
            hausdorff = d;
            result = j;
        }
    }
    return result;
}

VHandle IterLine::hausdorff(VHandle vh, VHandle spiral_vh) {
    VHandle result = spiral_vh;
    double hausdorff = (p(vh)-p(spiral_vh)).norm();

    for(VHandle j = next(spiral_vh); j.is_valid() ; j = next(j)) {
        double d = (p(vh) - p(j)).norm();
        if (d  < hausdorff) {
            hausdorff = d;
            result = j;
        }
    }
    return result;
}

double IterLine::length(IHandle ih) {
    VHandle s_vh = v_handle(ih);
    double l = 0;
    for(VHandle vh=next(s_vh);vh.i != s_vh.i;vh=next(vh)) {
        l += (p(vh) - p(prev(vh))).norm();
    }
    return l;
}

VHandle IterLine::before(VHandle vh, double length) {
    for(vh=prev(vh);;vh = prev(vh)) {
        length -= (p(vh) - p(next(vh))).norm();
        if(length < 0) {
            double alpha = - length / (p(vh) - p(next(vh))).norm();
            Ps.emplace_back(p(vh) * (1-alpha) + p(next(vh)) * alpha);
            Vs.push_back(VLinkNode{vh.i, next(vh).i, size_t(i_handle(vh).i)});
            link(next(vh)).prev = Vs.size()-1;
            link(vh).next = Vs.size()-1;
            return VHandle{int(Vs.size()-1)};
        }
    }
}

void IterLine::after_delete(VHandle s_vh, VHandle e_vh) {
    for(VHandle vh = s_vh; vh.i != e_vh.i;vh = next(vh)) {
        link(prev(vh)).next = -1;
        link(vh).prev = -1;
    }
    link(prev(e_vh)).next = -1;
    link(e_vh).prev = -1;
}

VHandle IterLine::after(VHandle vh, double length) {
    for(vh=next(vh);;) {
        length -= (p(vh) - p(prev(vh))).norm();
        if(length < 0) {
            double alpha = - length / (p(vh) - p(prev(vh))).norm();
            Ps.emplace_back(p(vh) * (1-alpha) + p(prev(vh)) * alpha);
            Vs.push_back(VLinkNode{ prev(vh).i, vh.i, size_t(i_handle(vh).i)});
            link(prev(vh)).next = Vs.size()-1;
            link(vh).prev = Vs.size()-1;
            return VHandle{int(Vs.size()-1)};
        }
        vh = next(vh);
    }
}

void IterLine::after_reverse(VHandle s_vh, VHandle e_vh) {
    std::swap(link(s_vh).prev, link(s_vh).next);
    for(VHandle vh = prev(s_vh); vh.i != e_vh.i; vh = prev(vh))
        std::swap(link(vh).prev, link(vh).next);
}

std::vector<IterLine::RegionNode> IterLine::region_extract(){
    std::vector<IterLine::RegionNode> regions;
    std::vector<int> i_flag(Is.size(), -1);
    for(int i = 0; i < Is.size(); i++) {
        if (Is[i].child != -1) continue;
        auto ih = IHandle{i};
        i_flag[ih.i] = regions.size();
        for(ih = parent(ih);ih.is_valid() and i_flag[ih.i] == -1; ih = parent(ih))
            i_flag[ih.i] = regions.size();

        regions.push_back({IHandle{i}, ih, ih.is_valid() ? i_flag[ih.i]: -1});
    }
    return regions;
}

bool IterLine::after_test(VHandle vh, VHandle next_vh) {
    for(vh = next(vh); vh.is_valid(); vh = next(vh))
        if(next_vh.i == vh.i) return true;
    return false;
}

std::tuple<VHandle,VHandle> IterLine::make_fermat_spiral(std::vector<RegionNode>& regions) {
    struct FermatSpiral{
        VHandle start, end;
        int parent_idx;
    };
    std::vector<FermatSpiral> fermat_spirals;
    for(auto& [s_ih, e_ih, parent_idx]: region_extract()) {
        auto [start, end] = spiral_from_maximum(s_ih, e_ih);
        fermat_spirals.push_back({start, end, parent_idx});
    }
    for(auto& fs: fermat_spirals) {
        if(fs.parent_idx == -1) continue;
        auto& parent = fermat_spirals[fs.parent_idx];
        VHandle p_start = hausdorff(fs.start, parent.start);
        VHandle p_end = hausdorff(fs.end, parent.start);

        if(!after_test(p_start, p_end)) {
            after_delete(next(p_end), p_start);
            after_reverse(fs.start, next(fs.end));

            link(p_start).prev = fs.start.i;
            link(fs.start).next = p_start.i;
            link(p_end).next = fs.end.i;
            link(fs.end).prev = p_end.i;

        }else {
            after_delete(next(p_start), p_end);
            link(p_start).next = fs.start.i;
            link(fs.start).prev = p_start.i;
            link(p_end).prev = fs.end.i;
            link(fs.end).next = p_end.i;
        }
    }

    return {fermat_spirals[0].start, fermat_spirals[0].end};
}

void IterLine::print_fermat_spiral(std::ostream& stream, VHandle vh){
    for(auto& p:Ps)
        stream << "v " << p[0] << ' ' << p[1]  <<' ' << p[2]  << std::endl;

    for(; next(vh).is_valid(); vh = next(vh))
        stream << "l " << vh.i+1 <<' ' << next(vh).i + 1 << std::endl;
}

std::tuple<VHandle,VHandle> IterLine::spiral_from_maximum(IHandle s_ih, IHandle e_ih){
    VHandle vh_out, after_vh;
    IHandle ih = s_ih;
    VHandle vh = v_handle(ih);
    bool should_reverse = true;
    while(parent(ih).i != e_ih.i) {
        VHandle downward_vh = hausdorff(vh, parent(ih));
        double dis = (p(downward_vh)-p(vh)).norm();

        if(!child(ih).is_valid()) {
            vh_out = before(vh, dis);
            after_delete(next(vh_out), vh);
        }

        after_vh = after(downward_vh, dis);
        auto before_vh = before(downward_vh, dis);
        auto upward_vh = before(before_vh, dis);
        after_delete(next(downward_vh), after_vh);
        after_delete(next(upward_vh), before_vh);

        if (should_reverse) {
            after_reverse(after_vh, next(upward_vh));

            link(downward_vh).next = vh.i;
            link(vh).prev = downward_vh.i;

            link(vh_out).next = upward_vh.i;
            link(upward_vh).prev = vh_out.i;

        }else {
            after_reverse(before_vh, next(downward_vh));

            link(downward_vh).prev = vh.i;
            link(vh).next = downward_vh.i;

            link(vh_out).prev = upward_vh.i;
            link(upward_vh).next = vh_out.i;
        }
        vh_out = before_vh;
        vh  = after_vh;
        ih = parent(ih);
        should_reverse = !should_reverse;
    }

    return should_reverse ? std::tuple{after_vh,vh_out}:std::tuple{vh_out,after_vh};
}