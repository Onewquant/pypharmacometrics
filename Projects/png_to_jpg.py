from PIL import Image
from pathlib import Path

# PNG 파일들이 있는 폴더
input_dir = Path(r"C:/Users/ilma0/Documents/AuraDoc/sources/에어비정관")

# JPG 저장 폴더 (없으면 자동 생성)
output_dir = input_dir / "jpg_converted"
output_dir.mkdir(exist_ok=True)

# 폴더 내 모든 png 순회
for png_file in input_dir.glob("*.png"):

    # 이미지 열기
    img = Image.open(png_file)

    # JPG는 투명도(alpha)를 지원하지 않으므로 RGB 변환
    rgb_img = img.convert("RGB")

    # 저장 경로
    jpg_path = output_dir / f"{png_file.stem}.jpg"

    # 고품질 저장
    rgb_img.save(
        jpg_path,
        format="JPEG",
        quality=100,      # 최대 품질
        subsampling=0,    # 색상 서브샘플링 최소화
        optimize=True
    )

    print(f"Saved: {jpg_path}")

print("완료")