#include <cstdint>
#include <vector>
#include <cstdlib>
#include <iostream>

template <typename T, std::size_t Alignment>
struct AlignedAllocator {
    using value_type = T;

    AlignedAllocator() noexcept = default;

    template <typename U>
    AlignedAllocator(const AlignedAllocator<U, Alignment>&) noexcept {}

    T* allocate(std::size_t n) {
        void* ptr = nullptr;
        if (posix_memalign(&ptr, Alignment, n * sizeof(T)) != 0) {
            throw std::bad_alloc();
        }
        return static_cast<T*>(ptr);
    }

    void deallocate(T* ptr, std::size_t) noexcept {
        free(ptr);
    }

    template <typename U>
    struct rebind {
        using other = AlignedAllocator<U, Alignment>;
    };
};

int main() {
    constexpr size_t alignment = 32;
    std::vector<float, AlignedAllocator<float, alignment>> vec(128);

    std::cout << "Vector address: " << static_cast<void*>(vec.data()) << std::endl;
    if (reinterpret_cast<uintptr_t>(vec.data()) % alignment == 0) {
        std::cout << "Memory is aligned!\n";
    } else {
        std::cout << "Memory is NOT aligned!\n";
    }

    return 0;
}
