//
// Created by plush-jill on 3/18/25.
//

#ifndef ALIGNED_ALLOCATOR_HPP
#define ALIGNED_ALLOCATOR_HPP

#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <format>

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

    static int check() {
        // constexpr size_t alignment = Alignment; // Например, для AVX (256 бит = 32 байта)
        std::vector<T, AlignedAllocator<T, Alignment>> vec(200);

        std::cout << "Vector address: " << static_cast<void*>(vec.data()) << std::endl;

        if (reinterpret_cast<uintptr_t>(vec.data()) % Alignment == 0) {
            std::cout << "Memory is aligned!\n";
        } else {
            std::cout << "Memory is NOT aligned!\n";
            return 1;
        }
        for (int i {}; i < 200; ++i) {
            std::cout << std::format("Memory is {} ({})\n",
                reinterpret_cast<uintptr_t>(vec.data() + i) % Alignment == 0 ? "aligned" : "not aligned",
                reinterpret_cast<uintptr_t>(vec.data() + i));
        }

        return 0;
    }
};



#endif //ALIGNED_ALLOCATOR_HPP
