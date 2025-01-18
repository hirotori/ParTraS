module geometry_m
    use kind_parameters_m
    implicit none
    
contains

    pure subroutine compute_cell_centers(cell2verts, verts, cell_centers)
        integer(ip),intent(in) :: cell2verts(:,:)
        real(dp),intent(in) :: verts(:,:)
        real(dp),intent(inout) :: cell_centers(:,:)

        integer ic, nv
        real(dp) x_(3,10)

        x_(:,:) = 0.d0
        do ic = 1, size(cell2verts, dim=2)
            nv = count(cell2verts(:,ic) > 0)
            x_(:,1:nv) = verts(:,cell2verts(1:nv,ic))
            cell_centers(:,ic) = sum(x_, dim=2)/dble(nv)
        end do

    end subroutine


    pure subroutine compute_face_centers_and_normals(face2verts, verts, face_centers, face_normals)
        integer(ip),intent(in) :: face2verts(:,:)
        real(dp),intent(in) :: verts(:,:)
        real(dp),intent(inout) :: face_centers(:,:)
        real(dp),intent(inout) :: face_normals(:,:)

        integer jf, nv
        real(dp) x_(3,10), x12_(3), x13_(3), n_(3)

        x_(:,:) = 0.d0
        do jf = 1, size(face2verts, dim=2)
            nv = count(face2verts(:,jf) > 0)
            x_(:,1:nv) = verts(:,face2verts(1:nv,jf))
            face_centers(:,jf) = sum(x_(:,1:nv), dim=2)/dble(nv)
            x12_(:) = x_(:,2) - x_(:,1)
            x13_(:) = x_(:,3) - x_(:,1)
            n_(:) = cross_prod_(x12_, x13_)
            face_normals(:,jf) = n_(:)/sqrt(sum(n_*n_))
        end do

    end subroutine

    pure function cross_prod_(a, b) result(cross)
        real(dp),intent(in) :: a(3)
        real(dp),intent(in) :: b(3)
        real(dp) cross(3)

        cross(1) = a(2)*b(3) - a(3)*b(2)
        cross(2) = a(3)*b(1) - a(1)*b(3)
        cross(3) = a(1)*b(2) - a(2)*b(1)

    end function


end module geometry_m